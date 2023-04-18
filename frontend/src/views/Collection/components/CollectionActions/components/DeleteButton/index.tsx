import { H6, Intent } from "@blueprintjs/core";
import loadable from "@loadable/component";
import { ReactElement, useState } from "react";
import { Collection } from "src/common/entities";
import { DeleteButton as Button } from "./style";
import { DeleteCollectionFn } from "src/views/Collection/components/CollectionActions";

const AsyncAlert = loadable(
  () =>
    /*webpackChunkName: 'src/components/Alert' */ import("src/components/Alert")
);

interface Props {
  collectionName: Collection["name"];
  disabled?: boolean;
  handleDeleteCollection: DeleteCollectionFn;
  isDeleting: boolean;
}

const DeleteCollectionButton = ({
  collectionName,
  disabled = false,
  handleDeleteCollection,
  isDeleting,
}: Props): ReactElement => {
  const [isOpen, setIsOpen] = useState(false);
  const handleHover = () => {
    AsyncAlert.preload();
  };

  const handleClick = () => {
    setIsOpen(!isOpen);
  };

  return (
    <>
      <Button
        data-testid="delete-collection-button"
        disabled={disabled}
        onClick={handleClick}
        onMouseEnter={handleHover}
        sdsStyle="square"
        sdsType="primary"
      >
        Delete
      </Button>
      {isOpen && (
        <AsyncAlert
          cancelButtonText={"Cancel"}
          confirmButtonText={`Delete Collection`}
          intent={Intent.DANGER}
          isOpen={isOpen}
          onCancel={handleClick}
          onConfirm={handleDeleteCollection}
          loading={isDeleting}
        >
          <>
            <H6>
              {`Are you sure you want to delete the ${collectionName}
                collection?`}
            </H6>
            <p>
              Datasets in this collection will no longer be available for
              download and associated CELLxGENE Discover visualizations will be
              deleted. You cannot undo this action.
            </p>
          </>
        </AsyncAlert>
      )}
    </>
  );
};

export default DeleteCollectionButton;
