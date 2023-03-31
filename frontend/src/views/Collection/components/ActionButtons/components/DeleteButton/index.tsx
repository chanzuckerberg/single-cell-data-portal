import { H6, Intent } from "@blueprintjs/core";
import loadable from "@loadable/component";
import { ReactElement, useState } from "react";
import { Collection } from "src/common/entities";
import { StyledOutlineButton } from "src/components/common/Button/common/style";

const AsyncAlert = loadable(
  () =>
    /*webpackChunkName: 'src/components/Alert' */ import("src/components/Alert")
);

interface Props {
  handleConfirm: () => void;
  collectionName: Collection["name"];
  loading: boolean;
}

const DeleteCollectionButton = ({
  handleConfirm,
  collectionName,
  loading,
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
      <StyledOutlineButton
        onClick={handleClick}
        text="Delete Collection"
        intent={Intent.DANGER}
        data-testid="delete-collection-button"
        minimal
        outlined
        onMouseEnter={handleHover}
      />
      {isOpen && (
        <AsyncAlert
          cancelButtonText={"Cancel"}
          confirmButtonText={`Delete Collection`}
          intent={Intent.DANGER}
          isOpen={isOpen}
          onCancel={handleClick}
          onConfirm={handleConfirm}
          loading={loading}
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
