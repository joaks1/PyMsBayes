#! /usr/bin/env python

import unittest
import os

from pymsbayes.utils import tempfs
from pymsbayes.utils.functions import mkdr
from pymsbayes.test.support import package_paths
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

class TempFileSystemTestCase(unittest.TestCase):
    def setUp(self):
        self.parent = package_paths.output_path()
        mkdr(self.parent)
        self.prefix = 'test-temp'

    def tearDown(self):
        os.rmdir(self.parent)

    def test_initiate_purge(self):
        tmp = tempfs.TempFileSystem(parent=self.parent, prefix=self.prefix)
        self.assertIsInstance(tmp, tempfs.TempFileSystem)
        self.assertTrue(os.path.exists(tmp.base_dir))
        self.assertTrue(os.path.isdir(tmp.base_dir))
        self.assertEqual(tmp.prefix, self.prefix)
        self.assertEqual(tmp.parent, self.parent.rstrip(os.path.sep))
        self.assertTrue(os.path.basename(tmp.base_dir).startswith(
                self.prefix))
        self.assertFalse(tmp.deleted)
        tmp.purge()
        self.assertTrue(tmp.deleted)

if __name__ == '__main__':
    unittest.main()

