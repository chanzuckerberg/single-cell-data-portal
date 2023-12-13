import { Dialog, IButtonProps } from "@blueprintjs/core";
import loadable from "@loadable/component";
import * as React from "react";
import { FC, useState } from "react";
import { QUERY_PARAMETERS } from "src/common/constants/queryParameters";
import { Collection } from "src/common/entities";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import { removeParams } from "src/common/utils/removeParams";
import { StyledButton } from "./style";
import { useRouter } from "next/router";

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

const CreateCollectionButton = (props: Partial<IButtonProps>) => (
  <StyledButton sdsStyle="square" sdsType="primary" {...props}>
    Create Collection
  </StyledButton>
);

const CreateCollection: FC<{
  className?: string;
  id?: Collection["id"];
  Button?: React.ElementType;
}> = ({ className, id, Button }) => {
  const router = useRouter();

  const isCurator = get(FEATURES.CURATOR) === BOOLEAN.TRUE;
  const urlParams = new URLSearchParams(window.location.search);
  const param = urlParams.get(QUERY_PARAMETERS.LOGIN_MODULE_REDIRECT);

  const shouldModuleOpen = param?.toLowerCase() === BOOLEAN.TRUE;

  const [isOpen, setIsOpen] = useState(shouldModuleOpen);
  const { data: userInfo, isLoading } = useUserInfo(isCurator);

  if (!isCurator || isLoading) {
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
  const OpenDialogButton = Button || CreateCollectionButton;

  return (
    <>
      <OpenDialogButton
        onMouseOver={() => config.content.preload()}
        onClick={toggleOpen}
        {...{ className }}
      />
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
      removeParams({ params: QUERY_PARAMETERS.LOGIN_MODULE_REDIRECT, router });
    }
  }
};

export default CreateCollection;
