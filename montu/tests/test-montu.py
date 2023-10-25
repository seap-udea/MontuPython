import unittest
from montu import *

class Test(unittest.TestCase):

    def test(self):
        print(np.pi)

    def test_load_kernels(self):  
        # Load all kernels
        Montu.load_kernels()

    def test_montime(self):
        # Input strings
        mtime = MonTime('-2500-01-01 12:00:00')
        mtime = MonTime('2501 b.c.e. 01-01 12:00:00.00')
        mtime = MonTime('bce2501-01-01 12:00:00.00')
        print(mtime)

        # Using mixed calendar
        mtime = MonTime('bce2501-01-22 12:00:00.00',calendar='mixed')
        print(mtime)

        # Test ranges
        mtime = MonTime('-1-01-01 12:00:00')
        mtime = MonTime('100-01-01 12:00:00')
        mtime = MonTime('1000-01-01 12:00:00')
        print(mtime)
        
        pass

if __name__=="__main__":
   unittest.main(argv=['first-arg-is-ignored'],exit=False)