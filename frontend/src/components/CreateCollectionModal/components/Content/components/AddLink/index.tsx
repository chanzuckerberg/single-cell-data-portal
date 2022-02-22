import { Menu, MenuItem, Popover, Position } from "@blueprintjs/core";
import * as React from "react";
import { FC } from "react";
import {
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
  COLLECTION_LINK_TYPE_OPTION_DOI_TEXT,
} from "src/common/entities";
import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";

type LinkOption = {
  text: string;
  value: COLLECTION_LINK_TYPE;
}; /* TODO(cc) relocate? */

interface Props {
  doiSelected: boolean;
  handleClick: (type: COLLECTION_LINK_TYPE) => void;
}

interface AddLinkProps extends Props {
  Button: React.ElementType;
  doiSelected: boolean;
  fill?: boolean;
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

/**
 * Returns options when filter feature flag is enabled:
 * - with updated display text for DOI link type to "Publication DOI", and
 * - omits DOI link type from the menu list if type is already selected.
 * TODO(cc) revert to filter function and update COLLECTION_LINK_TYPE_OPTIONS.DOI text to "Publication DOI" once filter feature flag is removed (#1718).
 * @param options
 * @param doiSelected
 */
const reduceOptions = (
  options: LinkOption[],
  doiSelected: boolean
): LinkOption[] => {
  return options.reduce((acc: LinkOption[], option): LinkOption[] => {
    if (option.value === COLLECTION_LINK_TYPE.DOI) {
      if (doiSelected) {
        return acc;
      } else {
        Object.assign(option, { text: COLLECTION_LINK_TYPE_OPTION_DOI_TEXT });
      }
    }
    acc.push(option);
    return acc;
  }, []);
};

const LinkTypes: FC<Props> = ({ doiSelected, handleClick }) => {
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const filteredOptions = isFilterEnabled
    ? reduceOptions(options, doiSelected)
    : options;
  return (
    <Menu>
      {filteredOptions.map(({ text, value }) => (
        <MenuItem key={value} text={text} onClick={() => handleClick(value)} />
      ))}
    </Menu>
  );
};

const AddLink: FC<AddLinkProps> = ({
  doiSelected,
  fill = false,
  handleClick,
  Button,
}) => {
  return (
    <Popover
      content={
        <LinkTypes doiSelected={doiSelected} handleClick={handleClick} />
      }
      fill={fill}
      position={Position.BOTTOM_LEFT}
      captureDismiss={true} // Not setting this led to this ridiculously obfuscated bug where choosing a link type would close the menu that held the edit collection button. Closing that menu would close the entire edit collection modal
    >
      <Button />
    </Popover>
  );
};

export default AddLink;
