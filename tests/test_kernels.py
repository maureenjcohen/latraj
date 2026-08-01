import custom_kernels

def test_kernels_importable():
    assert callable(custom_kernels.smagdiff)
    assert callable(custom_kernels.convection_ou)
