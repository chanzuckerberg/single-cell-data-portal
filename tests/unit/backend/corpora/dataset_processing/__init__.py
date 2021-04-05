import unittest
from .test_downloader import TestDownload
from .test_process import TestDatasetProcessing

dataset_processing_test_suite = unittest.TestSuite(tests=(TestDownload, TestDatasetProcessing))
