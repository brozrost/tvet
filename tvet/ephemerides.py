import requests
import csv
import io

HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"

class HorizonsError(RuntimeError):
    pass

def fetch_ephems(
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
        response = requests.get(HORIZONS_URL, params=params, timeout=timeout)
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

    return xyz

if __name__ == "__main__":
    o_xyz = fetch_ephems(
        body="499",
        center="500@399",
        start_time="2006-01-01",
        stop_time="2006-01-20",
        step_size="1 d",
    )

    s_xyz = fetch_ephems(
        body="499",
        center="500@10",
        start_time="2006-01-01",
        stop_time="2006-01-20",
        step_size="1 d",
    )

    print("o[0] XYZ (km):", o_xyz[0])
    print("s[0] XYZ (km):", s_xyz[0])
