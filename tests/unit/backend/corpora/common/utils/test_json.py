import datetime
import json
import unittest
from enum import Enum

from sqlalchemy import Column, String

from backend.corpora.common.corpora_orm import Base
from backend.corpora.common.entities.entity import Entity
from backend.corpora.common.utils.json import CustomJSONEncoder


class DBTest(Base):
    __tablename__ = "test"
    id = Column(String, primary_key=True)
    name = Column(String)


class TestCustomJSONEncoder(unittest.TestCase):
    def test_datetime(self):
        test_datetime_value = datetime.datetime.fromtimestamp(0)
        expected_datetime = str(test_datetime_value.timestamp())
        self._verify_json_encoding(test_datetime_value, expected_datetime)

    def test_timedelta(self):
        time_1 = datetime.datetime.fromtimestamp(0)
        time_2 = datetime.datetime.fromtimestamp(10)
        test_timedelta_value = time_1 - time_2
        expected_timedelta = f'"{test_timedelta_value}"'
        self._verify_json_encoding(test_timedelta_value, expected_timedelta)

    def test_enum(self):
        class EnumClass(Enum):
            TEST = "test"

        test_enum_value = EnumClass.TEST
        expected_enum = f'"{test_enum_value.name}"'
        self._verify_json_encoding(test_enum_value, expected_enum)

    def test_base(self):
        params = dict(id="foo", name=None)
        test_base = DBTest(**params)
        expected_base = json.dumps(params, sort_keys=True)
        self._verify_json_encoding(test_base, expected_base)
        self.assertDictEqual({k: v for k, v in test_base}, params)

    def test_entity(self):
        params = dict(id="foo", name="bar")
        test_entity = Entity(DBTest(**params))
        expected_entity = json.dumps(params, sort_keys=True)
        self._verify_json_encoding(test_entity, expected_entity)

    def test_unsupported_type(self):
        class Unsupported:
            foo = "bar"

        test_unsupported_type = Unsupported()
        with self.assertRaises(TypeError):
            json.dumps(test_unsupported_type, cls=CustomJSONEncoder)

    def test_illegal_characters(self):
        tests = [
            ("%", '"%"'),
            ("/", '"/"'),
            ('"', '"\\""'),
            ("'", '"\'"'),
            (
                "https://s3.us-west-2.amazonaws.com/bogus-bucket/test_generate_file_url.h5ad?AWSAccessKeyId"
                "=ASIA2F53T2CQNMN2LBN5&Signature=tohsICVhzs9sJdNrG06Yh8%2FpM8Y%3D&x-amz-security-token=FwoGZXIvYXdzEMf%"
                "2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaDI%2FAU5xheJUplKF8GiLmAj03%2BzMJ34xAA%2Bz3IQzXZEUdJQcw4kih0b74%2Fftiw%"
                "2BWBh6bFUsUdO60mZLWafKFatbegWbzxcuW7funOqkfXgL59s%2Bph8xnYOi1KDD88r%2FTiWsX9UZ1CgW9Fai%2FZYvqUKYJaRfAt"
                "Z9v6g9COoVSic5YITWjinfVv85XajO6IJvzwRTXKfXyVFGIrmVwarsUDyik2%2F0a3E9aSXLO9WPvM85%2FCCON00cZlubwab9O9R4"
                "ahZBTqJjvje4mWPfmqxAT1KEejSLEAE5Bv1OTUNaIYPYNM3nibLvIi%2B25YVd1WUgvdnxfWhTMMa%2FZH9qcr%2F2oir%2F%2BBIj"
                "pJUkVPnHLdh5%2BYdMaL1ABGq3%2BoMU4AbrmHZIc2U6uT7C3IoaZgk2KzWACZH7Li1DRoyiyurS2K3RPLUZ1eUKvwrGFWgs2ZNtpB"
                "aiQ39qepqVtAg1Ks9tc9pRJU3l1AVbaAGj1NMVtm6thR61w4AQNvwCoo27rF%2BgUyIz2Th4H%2B7PXBvEDNp8qoUYn7rMPp%2BFdE"
                "XQCO3V5kvabHZojO&Expires=1599772637",
                '"https://s3.us-west-2.amazonaws.com/bogus-bucket/test_generate_file_url.h5ad?AWSAccessKeyId=ASIA2F53T2'
                "CQNMN2LBN5&Signature=tohsICVhzs9sJdNrG06Yh8%2FpM8Y%3D&x-amz-security-token=FwoGZXIvYXdzEMf%2F%2F%2F%2"
                "F%2F%2F%2F%2F%2F%2FwEaDI%2FAU5xheJUplKF8GiLmAj03%2BzMJ34xAA%2Bz3IQzXZEUdJQcw4kih0b74%2Fftiw%2BWBh6bFUs"
                "UdO60mZLWafKFatbegWbzxcuW7funOqkfXgL59s%2Bph8xnYOi1KDD88r%2FTiWsX9UZ1CgW9Fai%2FZYvqUKYJaRfAtZ9v6g9COoV"
                "Sic5YITWjinfVv85XajO6IJvzwRTXKfXyVFGIrmVwarsUDyik2%2F0a3E9aSXLO9WPvM85%2FCCON00cZlubwab9O9R4ahZBTqJjvj"
                "e4mWPfmqxAT1KEejSLEAE5Bv1OTUNaIYPYNM3nibLvIi%2B25YVd1WUgvdnxfWhTMMa%2FZH9qcr%2F2oir%2F%2BBIjpJUkVPnHLd"
                "h5%2BYdMaL1ABGq3%2BoMU4AbrmHZIc2U6uT7C3IoaZgk2KzWACZH7Li1DRoyiyurS2K3RPLUZ1eUKvwrGFWgs2ZNtpBaiQ39qepqV"
                "tAg1Ks9tc9pRJU3l1AVbaAGj1NMVtm6thR61w4AQNvwCoo27rF%2BgUyIz2Th4H%2B7PXBvEDNp8qoUYn7rMPp%2BFdEXQCO3V5kva"
                'bHZojO&Expires=1599772637"',
            ),
        ]
        for test, expected in tests:
            with self.subTest(test):
                actual = json.dumps(test, cls=CustomJSONEncoder)
                self.assertEqual(expected, actual)

    def _verify_json_encoding(self, test_value, expected_value):
        actual_value = json.dumps(test_value, cls=CustomJSONEncoder, sort_keys=True)
        self.assertEqual(expected_value, actual_value)
