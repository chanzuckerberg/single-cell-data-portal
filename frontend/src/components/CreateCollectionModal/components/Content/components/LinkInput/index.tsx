import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { debounce } from "lodash";
import React, { FC, useState } from "react";
import {
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
} from "src/common/entities";
import Input from "src/components/common/Form/Input";
import { GRAY } from "src/components/common/theme";
import { DEBOUNCE_TIME_MS } from "../../common/constants";
import AddLink from "../AddLink";
import { IconWrapper, LinkWrapper, StyledButton } from "./style";

export type LinkValue = {
  id: number;
  value: string;
  isValid: boolean;
  index: number;
  name: string;
};

const DOI_PLACEHOLDER = "https://doi.org/10.1126/science.aax6234";
const LINK_PLACEHOLDER = "https://cellxgene.cziscience.com";

interface Props {
  handleChange: ({ id, value, isValid }: LinkValue) => void;
  id: number;
  index: number;
  type: COLLECTION_LINK_TYPE;
  handleDelete: (id: number) => void;
  defaultValue: string;
}

const LinkInput: FC<Props> = ({
  handleChange,
  handleDelete,
  id,
  type,
  index,
  defaultValue,
}) => {
  const option = COLLECTION_LINK_TYPE_OPTIONS[type];

  const { text, value } = option;

  const [name, setName] = useState("");
  const [linkType, setLinkType] = useState(value as COLLECTION_LINK_TYPE);

  const LinkTypeButton = () => (
    <Button minimal={true} rightIcon="caret-down">
      {COLLECTION_LINK_TYPE_OPTIONS[linkType].text}
    </Button>
  );

  return (
    <LinkWrapper>
      <AddLink handleClick={handleLinkTypeChange} Button={LinkTypeButton} />

      <Input
        name={value}
        text="Name"
        placeholder="Name"
        handleChange={handleNameChange}
      />
      <Input
        name={value}
        text="URL"
        syncValidation={[isValidHttpUrl, isDOILink(value)]}
        placeholder={
          value === COLLECTION_LINK_TYPE.DOI
            ? DOI_PLACEHOLDER
            : LINK_PLACEHOLDER
        }
        defaultValue={defaultValue}
        handleChange={debounce(handleChange_, DEBOUNCE_TIME_MS)}
      />
      <IconWrapper>
        <StyledButton
          minimal
          color={GRAY.A}
          icon={IconNames.CROSS}
          onClick={() => handleDelete(id)}
        />
      </IconWrapper>
    </LinkWrapper>
  );

  function handleNameChange({ value }: { isValid: boolean; value: string }) {
    setName(value);
  }
  function handleLinkTypeChange(newLinkType: COLLECTION_LINK_TYPE) {
    setLinkType(newLinkType);
  }

  function handleChange_({
    isValid: isValidFromInput,
    value,
  }: {
    isValid: boolean;
    value: string;
  }) {
    handleChange({
      id,
      index,
      isValid: isValidFromInput,
      name,
      value,
    });
  }
};

function isDOILink(
  type: COLLECTION_LINK_TYPE
): (value: string) => true | string {
  return (value: string) => {
    // Skip validation if type is not DOI
    if (type !== COLLECTION_LINK_TYPE.DOI) return true;

    const isValid = value.includes("doi.org");

    return isValid || "Please enter a valid DOI link";
  };
}

function isValidHttpUrl(url: string): true | string {
  let result;

  try {
    result = new URL(url);
  } catch {
    result = { protocol: "" };
  }

  const isValid = result.protocol === "http:" || result.protocol === "https:";

  return isValid || "Please enter a link";
}

export default LinkInput;
