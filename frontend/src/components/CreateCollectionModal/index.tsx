import { Dialog } from "@blueprintjs/core";
import React, { FC, useState } from "react";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import Content from "./components/Content";
import { StyledButton } from "./style";

const CreateCollection: FC = () => {
  const [isOpen, setIsOpen] = useState(false);

  if (get(FEATURES.CREATE_COLLECTION) !== BOOLEAN.TRUE) return null;

  const toggleOpen = () => setIsOpen(!isOpen);

  return (
    <>
      <StyledButton onClick={toggleOpen}>Create Collection</StyledButton>
      <Dialog
        isOpen={isOpen}
        onClose={toggleOpen}
        canEscapeKeyClose={false}
        canOutsideClickClose={false}
      >
        <Content onClose={toggleOpen} />
      </Dialog>
    </>
  );
};

export default CreateCollection;
