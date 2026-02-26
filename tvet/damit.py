import requests

DAMIT_URL = "https://damit.cuni.cz/generated_files/open/AsteroidModel"

class DamitError(RuntimeError):
    pass

def fetch_text(url: str, *, timeout: float) -> str:
    try:
        response = requests.get(url, timeout=timeout)
    except requests.RequestException as exc:
        raise DamitError(f"Network error: {exc}") from exc
    
    if response.status_code != 200:
        raise DamitError(f"HTTP {response.status_code}: {response.text}")
    
    return response.text

def fetch_obj(*, model_id: int | str, timeout: float = 30.0):
    url = f"{DAMIT_URL}/{model_id}/shape.obj"
    text = fetch_text(url, timeout=timeout)

    vertices = []
    faces = []

    for line in text.splitlines():
        if not line:
            continue

        if line.startswith("v"):
            _, x, y, z = line.split()
            vertices.append((float(x), float(y), float(z)))
        elif line.startswith("f"):
            _, i, j, k = line.split()
            faces.append((int(i) - 1, int(j) - 1, int(k) - 1))

    return vertices, faces

def fetch_spin(*, model_id: int | str, timeout: float = 30.0):
    url = f"{DAMIT_URL}/{model_id}/spin.txt"
    text = fetch_text(url, timeout=timeout)

    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
    if not lines:
        raise DamitError(f"spin.txt is empty or comment-only for model_id={model_id}")
    
    l, b, period = map(float, lines[0].split())
    epoch, phi0 = map(float, lines[1].split())
    scat = float(lines[2].split()[0])

    return l, b, period, epoch, phi0, scat

def main():
    model_id = 103
    vertices, faces = fetch_obj(model_id=model_id)
    spin = fetch_spin(model_id=model_id)

    print(len(vertices), len(faces))
    print(spin)

if __name__ == "__main__":
    main()