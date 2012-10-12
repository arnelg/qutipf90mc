#function to run the nose test scripts
def test():
    import nose
    nose.run(defaultTest="qutipf90mc.tests",argv=['nosetests', '-v']) #runs tests in qutipf90mc/tests module only

