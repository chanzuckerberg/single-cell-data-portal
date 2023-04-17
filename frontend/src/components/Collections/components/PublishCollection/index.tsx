import { H6, Intent } from "@blueprintjs/core";
import loadable from "@loadable/component";
import { Dispatch, FC, SetStateAction } from "react";
import { Collection } from "src/common/entities";
import Policy from "./components/Policy";
import { ActionButton as Button } from "src/views/Collection/components/CollectionActions/style";

const AsyncAlert = loadable(
  () =>
    /*webpackChunkName: 'src/components/Alert' */ import("src/components/Alert")
);

interface Props {
  handlePublishCollection: () => void;
  isPublishable: boolean;
  isPublishing: boolean;
  isPublishOpen: boolean;
  revision_of: Collection["revision_of"];
  setIsPublishOpen: Dispatch<SetStateAction<boolean>>;
}

const PublishCollection: FC<Props> = ({
  handlePublishCollection,
  isPublishable,
  isPublishing,
  isPublishOpen,
  revision_of,
  setIsPublishOpen,
}) => {
  const handleHover = () => {
    AsyncAlert.preload();
  };
  return (
    <>
      <Button
        data-testid="publish-collection-button"
        disabled={!isPublishable}
        onClick={() => setIsPublishOpen(true)}
        onMouseEnter={handleHover}
        sdsStyle="square"
        sdsType="primary"
      >
        Publish
      </Button>
      {isPublishOpen && (
        <AsyncAlert
          cancelButtonText={"Cancel"}
          confirmButtonText={
            revision_of ? "Publish Revision" : "Publish Collection"
          }
          intent={Intent.PRIMARY}
          isOpen={isPublishOpen}
          onCancel={() => setIsPublishOpen(false)}
          onConfirm={handlePublishCollection}
          loading={isPublishing}
        >
          {revision_of ? (
            <>
              <H6>
                Are you sure you want to publish a revision to this collection?
              </H6>
              <p>
                Any datasets you’ve removed will no longer be accessible to
                portal users. Links to deleted datasets from external sites will
                display a message that the dataset has been withdrawn by the
                publisher.
              </p>
            </>
          ) : (
            <>
              <H6>Are you sure you want to publish this collection?</H6>
              <p>
                The datasets, related metadata, and CXG(s) in this collection
                will be viewable and downloadable by all portal users. External
                sites may directly link to the CXG(s).
              </p>
            </>
          )}
          <Policy />
        </AsyncAlert>
      )}
    </>
  );
};

export default PublishCollection;
