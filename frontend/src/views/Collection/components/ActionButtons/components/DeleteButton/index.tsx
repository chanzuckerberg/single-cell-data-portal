import { Button, H6, Intent } from "@blueprintjs/core";
import loadable from "@loadable/component";
import { ReactElement, useState } from "react";
import { Collection } from "src/common/entities";
import { Wrapper } from "../../style";

interface Props {
  handleConfirm: () => void;
  collectionName: Collection["name"];
}

const DeleteCollectionButton = ({
  handleConfirm,
  collectionName,
}: Props): ReactElement => {
  const AsyncAlert = loadable(
    () =>
      /*webpackChunkName: 'src/components/Alert' */ import(
        "src/components/Alert"
      )
  );
  const [isOpen, setIsOpen] = useState(false);
  const handleHover = () => {
    AsyncAlert.preload();
  };

  const handleClick = () => {
    setIsOpen(!isOpen);
  };

  return (
    <Wrapper>
      <Button
        onClick={handleClick}
        text="Delete Collection"
        intent={Intent.DANGER}
        data-test-id="delete-collection-button"
        minimal
        outlined
        onMouseOver={handleHover}
      />
      {isOpen && (
        <AsyncAlert
          cancelButtonText={"Cancel"}
          confirmButtonText={`Delete Collection`}
          intent={Intent.DANGER}
          isOpen={isOpen}
          onCancel={handleClick}
          onConfirm={handleConfirm}
        >
          <>
            <H6>
              {`Are you sure you want to delete the ${collectionName}
                collection?`}
            </H6>
            <p>
              Datasets in this collection will no longer be available for
              download and associated cellxgene visualizations will be deleted.
              You cannot undo this action.
            </p>
          </>
        </AsyncAlert>
      )}
    </Wrapper>
  );
};

export default DeleteCollectionButton;
