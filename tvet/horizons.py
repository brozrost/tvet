import requests
import csv
import io
import numpy as np

HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"

class HorizonsError(RuntimeError):
    pass

class HorizonsClient:
    def __init__(self, *, base_url: str = HORIZONS_URL):
        self.base_url = base_url

    def _format_time(self, time: str) -> str:
        return time.replace("T", " ")

    def _send_request(self, *, params: dict, timeout: float) -> str:
        try:
            response = requests.get(self.base_url, params=params, timeout=timeout)
        except requests.RequestException as exc:
            raise HorizonsError(f"Network error: {exc}") from exc

        if response.status_code != 200:
            raise HorizonsError(f"HTTP {response.status_code}: {response.text}")

        try:
            data = response.json()
        except ValueError as exc:
            raise HorizonsError("Invalid JSON response from Horizons.") from exc

        if "error" in data and data["error"]:
            raise HorizonsError(data["error"])

        if "result" not in data:
            raise HorizonsError("Unexpected Horizons response format.")

        return data["result"]
    
    def _extract_lines(self, *, text: str) -> list[str]:
        start = text.find("$$SOE")
        end = text.find("$$EOE")
        if start == -1 or end == -1:
            raise HorizonsError("Ephemeris block not found in response.")

        block = text[start + 5:end].strip()
        lines = [ln for ln in block.splitlines() if ln.strip()]

        return lines
    
    def _parse_row(self, row: list[str]) -> tuple[float, float, float]:
        # Drop trailing comma
        if row and row[-1].strip() == "":
            row = row[:-1]

        if len(row) < 5:
            raise HorizonsError(f"Unexpected VECTORS row: {row}")

        try:
            x = float(row[2].strip())
            y = float(row[3].strip())
            z = float(row[4].strip())
        except ValueError as exc:
            raise HorizonsError(f"Failed to parse XYZ from row: {row}") from exc
        
        return x, y, z
    
    def _normalize_vectors(self, vectors: np.ndarray) -> np.ndarray:
        vectors = np.asarray(vectors, dtype=np.double)

        if vectors.ndim == 1 and vectors.shape[0] == 3:
            n = np.linalg.norm(vectors)
            return vectors / (n if n > 0.0 else 1.0)

        if vectors.ndim != 2 or vectors.shape[1] != 3:
            raise ValueError(f"Expected shape (N,3) or (3,), got {vectors.shape}")
        
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        norms = np.where(norms > 0.0, norms, 1.0)

        return vectors / norms

    def fetch_single_ephem(
        self,
        *,
        body: str,
        center: str,
        epoch: str,
        obj_data: bool = False,
        timeout: float = 30.0,
    ) -> np.ndarray:
        params = {
            "format": "json",
            "MAKE_EPHEM": "YES",
            "EPHEM_TYPE": "VECTORS",
            "COMMAND": f"'{body}'",
            "CENTER": f"'{center}'",
            "TLIST": f"'{self._format_time(epoch)}'",
            "VEC_TABLE": "'2'",
            "OUT_UNITS": "'KM-S'",
            "CSV_FORMAT": "YES",
            "OBJ_DATA": "'YES'" if obj_data else "'NO'",
        }

        response_text = self._send_request(params=params, timeout=timeout)

        lines = self._extract_lines(text=response_text)
        if len(lines) != 1:
            raise HorizonsError(f"Expected exactly one VECTORS row for epoch={epoch}, got {len(lines)}")
        
        reader = csv.reader(io.StringIO("\n".join(lines)))
        row = next(reader, None)
        if not row:
            raise HorizonsError("Empty VECTORS row in response.")
        
        x, y, z = self._parse_row(row)

        return np.array([x, y, z], dtype=np.double)

    def fetch_ephems(
        self,
        *,
        body: str,
        center: str,
        start_time: str,
        stop_time: str,
        step_size: str,
        obj_data: bool = False,
        timeout: float = 30.0,
    ) -> np.ndarray:
        params = {
            "format": "json",
            "MAKE_EPHEM": "YES",
            "EPHEM_TYPE": "VECTORS",
            "COMMAND": f"'{body}'",
            "CENTER": f"'{center}'",
            "START_TIME": f"'{start_time}'",
            "STOP_TIME": f"'{stop_time}'",
            "STEP_SIZE": f"'{step_size}'",
            "VEC_TABLE": "'2'",
            "OUT_UNITS": "'KM-S'",
            "CSV_FORMAT": "YES",
            "OBJ_DATA": "'YES'" if obj_data else "'NO'",
        }

        response_text = self._send_request(params=params, timeout=timeout)
        lines = self._extract_lines(text=response_text)

        reader = csv.reader(io.StringIO("\n".join(lines)))

        xyz = []
        for row in reader:
            if not row:
                continue

            x, y, z = self._parse_row(row)
            xyz.append((x, y, z))

        return np.asarray(xyz, dtype=np.double)
    
    def fetch_single_so(
        self,
        *,
        body: str,
        epoch: str,
        observer_center: str = "500@399",
        sun_center: str = "500@10",
        obj_data: bool = False,
        normalize: bool = True,
        timeout: float = 30.0,
    ):
        o = self.fetch_single_ephem(
            body=body,
            center=observer_center,
            epoch=epoch,
            obj_data=obj_data,
            timeout=timeout,
        )
    
        s = self.fetch_single_ephem(
            body=body,
            center=sun_center,
            epoch=epoch,
            obj_data=obj_data,
            timeout=timeout,
        )

        if o.shape != (3,) or s.shape != (3,):
            raise HorizonsError(f"Single vector shape mismatch: o={o.shape}, s={s.shape}")
        
        if not normalize:
            return s, o
        
        return self._normalize_vectors(s), self._normalize_vectors(o)

    def fetch_so(
        self,
        *,
        body: str,
        start_time: str,
        stop_time: str,
        step_size: str,
        observer_center: str = "500@399",
        sun_center: str = "500@10",
        obj_data: bool = False,
        normalize: bool = True,
        timeout: float = 30.0,
    ):
        o = self.fetch_ephems(
            body=body,
            center=observer_center,
            start_time=start_time,
            stop_time=stop_time,
            step_size=step_size,
            obj_data=obj_data,
            timeout=timeout,
        )

        s = self.fetch_ephems(
            body=body,
            center=sun_center,
            start_time=start_time,
            stop_time=stop_time,
            step_size=step_size,
            obj_data=obj_data,
            timeout=timeout,
        )

        if o.ndim != 2 or s.ndim != 2:
            raise HorizonsError(f"Ephemerides must be 2D arrays: o={o.shape}, s={s.shape}")
        if o.shape[1] != 3 or s.shape[1] != 3:
            raise HorizonsError(f"Ephemerides must have shape (N,3): o={o.shape}, s={s.shape}")
        if o.shape[0] == 0 or s.shape[0] == 0:
            raise HorizonsError("Empty ephemerides arrays returned.")
        if o.shape != s.shape:
            raise HorizonsError(f"Ephemerides shape mismatch: o={o.shape}, s={s.shape}")

        if not normalize:
            return s, o
        
        return self._normalize_vectors(s), self._normalize_vectors(o)
