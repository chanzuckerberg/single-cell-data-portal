import { Button, Intent, UL } from "@blueprintjs/core";
import React, { FC } from "react";
import EmptyModal from "../DatasetTab/components/EmptyModal";

const GenesetTab: FC = () => {
  return (
    <EmptyModal
      title="No gene sets uploaded"
      content={
        <div>
          Eget leo, urna sed at mattis nullam tincidunt.{" "}
          <UL>
            <li>Lorem ipsum dolor sit amet consectetur adipisicing elit.</li>
            <li>Facere libero doloremque nulla ex.</li>
          </UL>
        </div>
      }
      button={
        <Button intent={Intent.PRIMARY} minimal outlined text="Add Gene Set" />
      }
    />
  );
};

export default GenesetTab;
