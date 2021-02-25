from sqlalchemy import Column, String, ForeignKey
from sqlalchemy.orm import relationship

from backend.corpora.common.corpora_orm import Base
from backend.corpora.common.utils.db_session import DBSessionMaker, db_session_manager
from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class DbParent(Base):
    __tablename__ = "parent"

    name = Column(String)
    optional = Column(String)
    children = relationship("DbChild", back_populates="parent", cascade="all, delete-orphan")
    properties = relationship("DbProperty", secondary="associate_property_owner", back_populates="parents")
    #
    # many_to_many = relationship("Child", back_populates="genesets")
    # circular = relationship()


class DbChild(Base):
    __tablename__ = "child"

    name = Column(String)
    parent_id = Column(String, ForeignKey("parent.id"), index=True)
    parent = relationship("DbParent", uselist=False, back_populates="children")


class DbProperty(Base):
    __tablename__ = "property"
    name = Column(String)
    parents = relationship("DbParent", secondary="associate_property_owner", back_populates="properties")


class AssociatePropertyOwner(Base):
    __tablename__ = "associate_property_owner"
    parent_id = Column(String, ForeignKey("parent.id"), index=True)
    property_id = Column(String, ForeignKey("property.id"), index=True)


# parent:1->N:child
# parent:N->M:property
# child:1->M:property


class testToDict(DataPortalTestCase):
    @classmethod
    def setUpClass(cls):
        DataPortalTestCase.setUpClass()
        engine = DBSessionMaker().engine
        Base.metadata.drop_all(engine)
        Base.metadata.create_all(engine)

    @classmethod
    def tearDownClass(cls):
        cls.reinitialize_database()
        DataPortalTestCase.tearDownClass()

    @staticmethod
    def remove_ids(new_dict) -> dict:
        def _remove_ids(_new_dict):
            if not isinstance(_new_dict, dict):
                return
            _new_dict.pop("id", None)
            for value in _new_dict.values():
                if isinstance(value, dict):
                    _remove_ids(value)
                elif isinstance(value, list):
                    for list_value in value:
                        _remove_ids(list_value)
            return _new_dict

        return _remove_ids(new_dict)

    @staticmethod
    def db_delete(table, _id):
        with db_session_manager() as session:
            col = session.query(table).get(_id)
            if col:
                session.delete(col)

    def create_db_object(self, table, **params):
        db_obj = table(**params)
        self.session.add(db_obj)
        self.session.commit()
        # Cleanup collection after test
        self.addCleanup(self.db_delete, table, db_obj.id)
        return db_obj

    def test_defaults(self):
        parent = self.create_db_object(DbParent, name="bob")
        result = parent.to_dict()
        result = self.remove_ids(result)
        self.assertDictEqual({"name": "bob", "children": [], "properties": [], "optional": None}, result)

    def test_remove_none(self):
        parent = self.create_db_object(DbParent, name="bob")
        result = parent.to_dict(remove_none=True)
        result = self.remove_ids(result)
        self.assertDictEqual({"name": "bob", "children": [], "properties": []}, result)

    def test_children(self):
        parent = self.create_db_object(DbParent, name="bob", children=[DbChild(name="alice")])
        result = parent.to_dict(remove_none=True)
        result = self.remove_ids(result)
        self.assertDictEqual({"name": "bob", "children": [{"name": "alice", "parent_id": parent.id}], "properties": []},
                             result)

    def test_remove_attr(self):
        parent = self.create_db_object(DbParent, name="bob")
        result = parent.to_dict(remove_attr=["children", "properties", "optional"])
        result = self.remove_ids(result)
        self.assertDictEqual({"name": "bob"},
                             result)
