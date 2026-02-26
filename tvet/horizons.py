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
    ):
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

        result_text = data["result"]

        start = result_text.find("$$SOE")
        end = result_text.find("$$EOE")
        if start == -1 or end == -1:
            raise HorizonsError("Ephemeris block not found in response.")

        block = result_text[start + 5:end].strip()

        lines = [ln for ln in block.splitlines() if ln.strip()]
        reader = csv.reader(io.StringIO("\n".join(lines)))

        xyz = []
        for row in reader:
            if not row:
                continue

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

            xyz.append((x, y, z))

        return np.asarray(xyz, dtype=np.double)
    
    def normalize_vectors(self, vectors: np.ndarray) -> np.ndarray:
        vectors = np.asarray(vectors, dtype=np.double)

        if vectors.ndim == 1 and vectors.shape[0] == 3:
            n = np.linalg.norm(vectors)
            return vectors / (n if n > 0.0 else 1.0)

        if vectors.ndim != 2 or vectors.shape[1] != 3:
            raise ValueError(f"Expected shape (N,3) or (3,), got {vectors.shape}")
        
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        norms = np.where(norms > 0.0, norms, 1.0)

        return vectors / norms
    
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

        if o.shape != s.shape:
            raise HorizonsError(f"Ephemerides shape mismatch: o={o.shape}, s={s.shape}")

        if not normalize:
            return s, o
        
        return self.normalize_vectors(s), self.normalize_vectors(o)
