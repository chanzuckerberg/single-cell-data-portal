import React, { FC, useEffect } from "react";
import dropins from "./dropins";

// (thuang): My personal Dropbox app key. Will be replaced with
// Single Cell app key once the contract goes through
const TEMP_APP_KEY = "66mv8tdrtxi2s0m";

// (thuang): `dropins.ts` library looks for this element id
const ID = "dropboxjs";

const Chooser: FC = () => {
  useEffect(() => {
    dropins();
  }, []);

  return <div id={ID} data-app-key={TEMP_APP_KEY}></div>;
};

export default Chooser;
