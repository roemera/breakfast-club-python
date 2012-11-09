import unittest
from barnard import Barnard

class BarnardTests(unittest.TestCase):
    def setUp(self):
        self.barnard = Barnard()
    def test(self):
        results = self.barnard.test(7,12,8,3)
        self.failUnless(results == 0.0341)

def main():
    unittest.main()

if __name__ == '__main__':
    main()