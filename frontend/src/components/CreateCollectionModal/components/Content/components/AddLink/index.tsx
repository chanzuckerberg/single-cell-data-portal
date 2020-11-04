import {
  Button,
  Intent,
  Menu,
  MenuItem,
  Popover,
  Position,
} from "@blueprintjs/core";
import React, { FC } from "react";
import { OPTIONS, TYPES } from "../LinkInput";

interface Props {
  handleClick: (type: TYPES) => void;
}

const OPTION_ORDER = [
  TYPES.DOI,
  TYPES.RAW_DATA,
  TYPES.PROTOCOL,
  TYPES.LAB_WEBSITE,
  TYPES.OTHER,
];

const options = OPTION_ORDER.map((type) => OPTIONS[type]);

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
      <Button minimal intent={Intent.PRIMARY}>
        Add Link
      </Button>
    </Popover>
  );
};

export default AddLink;
