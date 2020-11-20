python_env <- BasiliskEnvironment(
    "densvis",
    pkgname = "densvis",
    channels = c("conda-forge", "alanocallaghan"),
    packages = c(
        "densmap-learn=0.2.2",
        "scikit-learn=0.23.2",
        "numba=0.48.0",
        "numpy=1.19.1"
    )
)
