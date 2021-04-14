import { Dialog } from "@blueprintjs/core";
import loadable from "@loadable/component";
import React, { FC, useState } from "react";
import { QUERY_PARAMETERS } from "src/common/constants/queryParameters";
import { Collection } from "src/common/entities";
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

const CreateCollection: FC<{ className?: string; id?: Collection["id"] }> = ({
  className,
  id,
}) => {
  const isAuth = get(FEATURES.AUTH) === BOOLEAN.TRUE;
  const urlParams = new URLSearchParams(window.location.search);
  const param = urlParams.get(QUERY_PARAMETERS.LOGIN_MODULE_REDIRECT);

  const shouldModuleOpen = param?.toLowerCase() === BOOLEAN.TRUE;

  const [isOpen, setIsOpen] = useState(shouldModuleOpen);
  const { data: userInfo, isLoading } = useUserInfo(isAuth);

  if (get(FEATURES.CREATE_COLLECTION) !== BOOLEAN.TRUE || isLoading) {
    return null;
  }

  const config = userInfo?.name
    ? {
        canEscapeKeyClose: false,
        canOutsideClickClose: false,
        content: AsyncContent,
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
        {...{ className }}
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
        {isOpen && <config.content onClose={toggleOpen} id={id} />}
      </Dialog>
    </>
  );

  function toggleOpen() {
    setIsOpen(!isOpen);
    if (shouldModuleOpen) {
      const url = window.location.href;
      const afterSlashBeforeParam = url
        .substring(url.indexOf("/") + 1)
        .split("?")[0];

      urlParams.delete(QUERY_PARAMETERS.LOGIN_MODULE_REDIRECT);
      if (urlParams.toString().length > 0) {
        const newURL = afterSlashBeforeParam + "?" + urlParams.toString();
        window.history.replaceState(null, " ", "/" + newURL);
      }
    }
  }
};

export default CreateCollection;
