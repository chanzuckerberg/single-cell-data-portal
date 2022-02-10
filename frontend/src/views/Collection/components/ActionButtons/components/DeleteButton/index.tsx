import { Button, H6, Intent } from "@blueprintjs/core";
import loadable from "@loadable/component";
import { Fragment, ReactElement, useState } from "react";
import { Collection } from "src/common/entities";
import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { Wrapper } from "../../style";

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
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const Actions = isFilterEnabled ? Fragment : Wrapper;
  const [isOpen, setIsOpen] = useState(false);
  const handleHover = () => {
    AsyncAlert.preload();
  };

  const handleClick = () => {
    setIsOpen(!isOpen);
  };

  return (
    <Actions>
      <Button
        onClick={handleClick}
        text="Delete Collection"
        intent={Intent.DANGER}
        data-test-id="delete-collection-button"
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
              download and associated cellxgene visualizations will be deleted.
              You cannot undo this action.
            </p>
          </>
        </AsyncAlert>
      )}
    </Actions>
  );
};

export default DeleteCollectionButton;
