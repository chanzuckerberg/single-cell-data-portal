import { IconNames } from "@blueprintjs/icons";
import { debounce } from "lodash";
import React, { FC } from "react";
import {
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
} from "src/common/entities";
import Input from "src/components/common/Form/Input";
import { LabelText, StyledDiv } from "src/components/common/Form/Input/style";
import { GRAY } from "src/components/common/theme";
import { DEBOUNCE_TIME_MS } from "../../common/constants";
import AddLink from "../AddLink";
import {
  IconWrapper,
  LinkWrapper,
  StyledButton,
  StyledLinkTypeButton,
  StyledURLInput,
} from "./style";

export type LinkValue = {
  id: number;
  url: string;
  isValid: boolean;
  index: number;
  linkName: string;
  linkType: COLLECTION_LINK_TYPE;
};

const DOI_PLACEHOLDER = "https://doi.org/10.1126/science.aax6234";
const LINK_PLACEHOLDER = "https://cellxgene.cziscience.com";

interface Props {
  handleChange: ({ id, url, isValid }: LinkValue) => void;
  id: number;
  index: number;
  linkName: string;
  linkType: COLLECTION_LINK_TYPE;
  handleDelete: (id: number) => void;
  url: string;
  isValid: boolean;
}

const LinkInput: FC<Props> = ({
  handleChange,
  handleDelete,
  id,
  linkType,
  index,
  url,
  linkName,
  isValid,
}) => {
  const option = COLLECTION_LINK_TYPE_OPTIONS[linkType];
  const { text, value } = option;

  const LinkTypeButton = () => (
    <StyledDiv>
      <LabelText>Type</LabelText>
      <StyledLinkTypeButton outlined minimal={true} rightIcon="caret-down">
        {text}
      </StyledLinkTypeButton>
    </StyledDiv>
  );

  return (
    <LinkWrapper>
      <AddLink handleClick={handleLinkTypeChange} Button={LinkTypeButton} />
      <Input
        name="Name"
        text="Name(optional)"
        placeholder="Name"
        handleChange={handleNameChange}
        defaultValue={linkName}
        percentage={25}
      />
      <StyledURLInput
        name={value}
        text="URL"
        syncValidation={[isValidHttpUrl, isDOILink(value)]}
        placeholder={
          value === COLLECTION_LINK_TYPE.DOI
            ? DOI_PLACEHOLDER
            : LINK_PLACEHOLDER
        }
        defaultValue={url}
        handleChange={debounce(handleChange_, DEBOUNCE_TIME_MS)}
        percentage={40}
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

  function handleNameChange({
    value: newLinkName,
  }: {
    isValid: boolean;
    value: string;
  }) {
    handleChange({
      id,
      index,
      isValid,
      linkName: newLinkName,
      linkType,
      url,
    });
  }
  function handleLinkTypeChange(newLinkType: COLLECTION_LINK_TYPE) {
    handleChange({
      id,
      index,
      isValid,
      linkName,
      linkType: newLinkType,
      url,
    });
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
      linkName,
      linkType,
      url: value,
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
