from . import conversions

# MARK: - tlist_format()

def tlist_format(value: str) -> str:
    value = value.strip()

    if not value:
        raise ValueError("Epoch string is empty.")

    if value[:2].upper() == "JD":
        jd = float(value[2:])

        return conversions.jd_to_iso(jd).replace("T", " ")

    try:
        jd = float(value)

        return conversions.jd_to_iso(jd).replace("T", " ")
    
    except ValueError:
        pass

    return value.replace("T", " ")