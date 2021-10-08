import { Button, Intent, UL } from "@blueprintjs/core";
import { FC } from "react";
import { StyledLink } from "../../common/style";
import EmptyModal from "../EmptyModal";

// @seve TODO: Link needs to be replaced when we have the format
const CLI_README_LINK =
  "https://github.com/chanzuckerberg/single-cell-curation/blob/main/readme.md";

const GeneSetTab: FC = () => {
  return (
    <EmptyModal
      title="No gene sets uploaded"
      content={
        <div>
          When uploading gene sets, please keep in mind
          <UL>
            <li>
              You must upload gene sets using{" "}
              <StyledLink href={CLI_README_LINK}>this format.</StyledLink>
            </li>
            <li>
              You may associate gene sets with datasets in the dataset tab.
            </li>
          </UL>
        </div>
      }
      button={
        <Button intent={Intent.PRIMARY} minimal outlined text="Add Gene Set" />
      }
    />
  );
};

export default GeneSetTab;
