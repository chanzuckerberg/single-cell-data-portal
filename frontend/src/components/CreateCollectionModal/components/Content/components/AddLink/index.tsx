import {
  Button,
  Intent,
  Menu,
  MenuItem,
  Popover,
  Position,
} from "@blueprintjs/core";
import React, { FC } from "react";
import {
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
} from "src/common/entities";

interface Props {
  handleClick: (type: COLLECTION_LINK_TYPE) => void;
}

const OPTION_ORDER = [
  COLLECTION_LINK_TYPE.DOI,
  COLLECTION_LINK_TYPE.RAW_DATA,
  COLLECTION_LINK_TYPE.PROTOCOL,
  COLLECTION_LINK_TYPE.LAB_WEBSITE,
  COLLECTION_LINK_TYPE.OTHER,
];

const options = OPTION_ORDER.map((type) => COLLECTION_LINK_TYPE_OPTIONS[type]);

const LinkTypes: FC<Props> = ({ handleClick }) => {
  return (
    <Menu>
      {options.map(({ text, value }) => (
        <MenuItem key={value} text={text} onClick={() => handleClick(value)} />
      ))}
    </Menu>
  );
};

const AddLink: FC<Props> = ({ handleClick }) => {
  return (
    <Popover
      content={<LinkTypes handleClick={handleClick} />}
      position={Position.BOTTOM_LEFT}
    >
      <Button outlined intent={Intent.PRIMARY}>
        Add Link
      </Button>
    </Popover>
  );
};

export default AddLink;
