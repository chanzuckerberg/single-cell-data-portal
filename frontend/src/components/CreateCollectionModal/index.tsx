import { Dialog } from "@blueprintjs/core";
import loadable from "@loadable/component";
import React, { FC, useState } from "react";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import { StyledButton } from "./style";

const AsyncContent = loadable(
  () =>
    /*webpackChunkName: 'CreateCollectionModalContent' */ import(
      "./components/Content"
    )
);

const AsyncCTA = loadable(
  () =>
    /*webpackChunkName: 'CreateCollectionModalCTA' */ import("./components/CTA")
);

const CreateCollection: FC = () => {
  const isAuth = get(FEATURES.AUTH) === BOOLEAN.TRUE;

  const [isOpen, setIsOpen] = useState(false);
  const { data: userInfo, isLoading } = useUserInfo(isAuth);

  if (get(FEATURES.CREATE_COLLECTION) !== BOOLEAN.TRUE || isLoading) {
    return null;
  }

  const config = userInfo?.name
    ? {
        canEscapeKeyClose: false,
        canOutsideClickClose: false,
        content: AsyncContent,
        isCloseButtonShown: false,
      }
    : {
        canEscapeKeyClose: true,
        canOutsideClickClose: true,
        content: AsyncCTA,
        isCloseButtonShown: true,
        title: "Create an account or sign-in to get started",
      };

  return (
    <>
      <StyledButton
        onMouseOver={() => config.content.preload()}
        onClick={toggleOpen}
      >
        Create Collection
      </StyledButton>
      <Dialog
        isCloseButtonShown={config.isCloseButtonShown}
        title={config.title}
        isOpen={isOpen}
        onClose={toggleOpen}
        canEscapeKeyClose={config.canEscapeKeyClose}
        canOutsideClickClose={config.canOutsideClickClose}
      >
        {isOpen && <config.content onClose={toggleOpen} />}
      </Dialog>
    </>
  );

  function toggleOpen() {
    setIsOpen(!isOpen);
  }
};

export default CreateCollection;
