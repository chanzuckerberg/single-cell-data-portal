import React from "react";
import ReactDOM from "react-dom";
import { Provider } from "react-redux";
import { HotkeysProvider, FocusStyleManager } from "@blueprintjs/core";

import "./index.css";

/* our code */
import App from "./components/app";
import store from "./reducers";

FocusStyleManager.onlyShowFocusOnTabs();

ReactDOM.render(
  <Provider store={store}>
    <HotkeysProvider>
      <App />
    </HotkeysProvider>
  </Provider>,
  document.getElementById("root"),
);
