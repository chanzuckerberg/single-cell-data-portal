import { Menu, MenuItem, Popover, Position } from "@blueprintjs/core";
import * as React from "react";
import { FC } from "react";
import {
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
} from "src/common/entities";

interface Props {
  handleClick: (type: COLLECTION_LINK_TYPE) => void;
}

interface AddLinkProps extends Props {
  Button: React.ElementType;
}

const OPTION_ORDER = [
  COLLECTION_LINK_TYPE.DOI,
  COLLECTION_LINK_TYPE.DATA_SOURCE,
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

const AddLink: FC<AddLinkProps> = ({ handleClick, Button }) => {
  return (
    <Popover
      content={<LinkTypes handleClick={handleClick} />}
      position={Position.BOTTOM_LEFT}
      captureDismiss={true} // Not setting this led to this ridiculously obfuscated bug where choosing a link type would close the menu that held the edit collection button. Closing that menu would close the entire edit collection modal
    >
      <Button />
    </Popover>
  );
};

export default AddLink;
