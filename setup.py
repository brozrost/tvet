from setuptools import Extension, setup # type: ignore

setup(
    ext_modules=[
        Extension(
            "tvet._tvet",
            sources=[
                "tvet/_tvet.c",
                "tvet/csrc/vector_math.c",
                "tvet/csrc/intersect_AB_t.c",
                "tvet/csrc/bounding_box.c",
                "tvet/csrc/shadowing.c",
            ],
            include_dirs=["tvet/csrc"],
        )
    ]
)