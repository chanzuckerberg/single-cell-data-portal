import { Dialog } from "@blueprintjs/core";
import loadable from "@loadable/component";
import React, { FC, useState } from "react";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { StyledButton } from "./style";

const AsyncContent = loadable(
  () =>
    /*webpackChunkName: 'CreateCollectionModalContent' */ import(
      "./components/Content"
    )
);

const CreateCollection: FC = () => {
  const [isOpen, setIsOpen] = useState(false);

  if (get(FEATURES.CREATE_COLLECTION) !== BOOLEAN.TRUE) return null;

  const toggleOpen = () => setIsOpen(!isOpen);

  return (
    <>
      <StyledButton
        onMouseOver={() => AsyncContent.preload()}
        onClick={toggleOpen}
      >
        Create Collection
      </StyledButton>
      <Dialog
        isOpen={isOpen}
        onClose={toggleOpen}
        canEscapeKeyClose={false}
        canOutsideClickClose={false}
      >
        {isOpen && <AsyncContent onClose={toggleOpen} />}
      </Dialog>
    </>
  );
};

export default CreateCollection;
