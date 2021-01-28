import React, { FC, useEffect } from "react";
import dropins from "./dropins";

// (thuang): Dropbox Business doesn't have shared account. So I'm using
// my work Dropbox's app key
const APP_KEY = "dernv4sjibm74ck";

// (thuang): `dropins.ts` library looks for this element id
const ID = "dropboxjs";

const Chooser: FC = () => {
  useEffect(() => {
    dropins();
  }, []);

  return <div id={ID} data-app-key={APP_KEY}></div>;
};

export default Chooser;
