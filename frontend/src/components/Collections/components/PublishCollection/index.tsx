import { Button, H6, Intent } from "@blueprintjs/core";
import loadable from "@loadable/component";
import { useRouter } from "next/router";
import React, { FC, useState } from "react";
import { ROUTES } from "src/common/constants/routes";
import { Collection } from "src/common/entities";
import { usePublishCollection } from "src/common/queries/collections";

const AsyncAlert = loadable(
  () =>
    /*webpackChunkName: 'src/components/Alert' */ import("src/components/Alert")
);

interface Props {
  id: Collection["id"];
  isPublishable: boolean;
}

const PublishCollection: FC<Props> = ({ id = "", isPublishable }) => {
  const [isOpen, setIsOpen] = useState(false);
  const [publish, { isSuccess, isLoading }] = usePublishCollection();
  const router = useRouter();

  if (isSuccess) {
    router.push(ROUTES.COLLECTION.replace(":id", id));
  }

  const toggleAlert = () => setIsOpen(!isOpen);

  const handleConfirm = async () => {
    await publish(id);
    toggleAlert();
  };

  const handleHover = () => {
    AsyncAlert.preload();
  };

  const handleClick = () => {
    setIsOpen(!isOpen);
  };

  return (
    <>
      <Button
        onMouseEnter={handleHover}
        onClick={handleClick}
        intent={Intent.PRIMARY}
        text="Publish"
        disabled={!isPublishable}
        data-test-id="publish-collection-button"
      />
      {isOpen && (
        <AsyncAlert
          cancelButtonText={"Cancel"}
          confirmButtonText={"Publish Collection"}
          intent={Intent.PRIMARY}
          isOpen={isOpen}
          onCancel={toggleAlert}
          onConfirm={handleConfirm}
          loading={isLoading}
        >
          <H6>Are you sure you want to publish this collection?</H6>
          <p>
            The datasets, related metadata, and CXG(s) in this collection will
            be viewable and downloadable by all portal users. External sites may
            directly link to the CXG(s).
          </p>
          <p>This action cannot be undone without manual intervention.</p>
        </AsyncAlert>
      )}
    </>
  );
};

export default PublishCollection;
