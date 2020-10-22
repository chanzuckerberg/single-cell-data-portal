import loadable from "@loadable/component";
import { Router } from "@reach/router";
import React, { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";

const AsyncHomepage = loadable(
  () => /*webpackChunkName: 'Homepage' */ import("src/views/Homepage")
);

const AsyncMyCollections = loadable(
  () =>
    /*webpackChunkName: 'AsyncMyCollections' */ import(
      "src/views/MyCollections"
    )
);

const AsyncCollection = loadable(
  () => /*webpackChunkName: 'AsyncCollection' */ import("src/views/Collection")
);

// (thuang): Check for Create Collection feature flag
get(FEATURES.CREATE_COLLECTION);

const AppContainer: FC = () => {
  return (
    <Router>
      <AsyncHomepage path={ROUTES.HOMEPAGE} />
      <AsyncMyCollections path={ROUTES.MY_COLLECTIONS} />
      <AsyncCollection path={ROUTES.COLLECTION} />
    </Router>
  );
};

export default AppContainer;
