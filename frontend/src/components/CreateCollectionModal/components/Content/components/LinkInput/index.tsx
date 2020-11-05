import { IconNames } from "@blueprintjs/icons";
import { debounce } from "lodash-es";
import React, { FC } from "react";
import {
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
} from "src/common/entities";
import Input from "src/components/common/Form/Input";
import { GRAY } from "src/components/common/theme";
import { DEBOUNCE_TIME_MS } from "src/components/CreateCollectionModal/components/Content/common/constants";
import { IconWrapper, StyledButton, Wrapper } from "./style";

export type LinkValue = {
  id: number;
  value: string;
  isValid: boolean;
  index: number;
};

const DOI_PLACEHOLDER = "https://doi.org/10.1126/science.aax6234";
const LINK_PLACEHOLDER = "https://cellxgene.cziscience.com";

interface Props {
  handleChange: ({ id, value, isValid }: LinkValue) => void;
  id: number;
  index: number;
  type: COLLECTION_LINK_TYPE;
  handleDelete: (id: number) => void;
}

const LinkInput: FC<Props> = ({
  handleChange,
  handleDelete,
  id,
  type,
  index,
}) => {
  const option = COLLECTION_LINK_TYPE_OPTIONS[type];

  const { text, value } = option;

  return (
    <Wrapper>
      <Input
        name={value}
        text={text}
        handleChange={debounce(handleChange_, DEBOUNCE_TIME_MS)}
        syncValidation={[isValidHttpUrl, isDOILink(value)]}
        noNameAttr
        placeholder={
          value === COLLECTION_LINK_TYPE.DOI
            ? DOI_PLACEHOLDER
            : LINK_PLACEHOLDER
        }
      />
      <IconWrapper>
        <StyledButton
          minimal
          color={GRAY.A}
          icon={IconNames.CROSS}
          onClick={() => handleDelete(id)}
        />
      </IconWrapper>
    </Wrapper>
  );

  function handleChange_({
    isValid: isValidFromInput,
    value,
  }: {
    isValid: boolean;
    value: string;
  }) {
    handleChange({ id, index, isValid: isValidFromInput, value });
  }
};

function isDOILink(
  type: COLLECTION_LINK_TYPE
): (value: string) => true | string {
  return (value: string) => {
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
