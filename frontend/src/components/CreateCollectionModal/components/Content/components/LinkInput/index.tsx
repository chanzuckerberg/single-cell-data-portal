import { IconNames } from "@blueprintjs/icons";
import { debounce } from "lodash-es";
import React, { FC } from "react";
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

export enum TYPES {
  DOI = "DOI",
  RAW_DATA = "RAW_DATA",
  PROTOCOL = "PROTOCOL",
  LAB_WEBSITE = "LAB_WEBSITE",
  OTHER = "OTHER",
}

export const OPTIONS = {
  [TYPES.DOI]: { text: "DOI", value: TYPES.DOI },
  [TYPES.RAW_DATA]: { text: "Raw Data", value: TYPES.RAW_DATA },
  [TYPES.PROTOCOL]: { text: "Protocol", value: TYPES.PROTOCOL },
  [TYPES.LAB_WEBSITE]: { text: "Lab Website", value: TYPES.LAB_WEBSITE },
  [TYPES.OTHER]: { text: "Other", value: TYPES.OTHER },
};

const DOI_PLACEHOLDER = "https://doi.org/10.1126/science.aax6234";
const LINK_PLACEHOLDER = "https://cellxgene.cziscience.com";

interface Props {
  handleChange: ({ id, value, isValid }: LinkValue) => void;
  id: number;
  index: number;
  type: TYPES;
  handleDelete: (id: number) => void;
}

const LinkInput: FC<Props> = ({
  handleChange,
  handleDelete,
  id,
  type,
  index,
}) => {
  const option = OPTIONS[type];

  const { text, value } = option;

  return (
    <Wrapper>
      <Input
        name={value}
        text={text}
        handleChange={debounce(handleChange_, DEBOUNCE_TIME_MS)}
        syncValidation={[isValidHttpUrl, isDOILink(value)]}
        noNameAttr
        placeholder={value === TYPES.DOI ? DOI_PLACEHOLDER : LINK_PLACEHOLDER}
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

function isDOILink(type: TYPES): (value: string) => true | string {
  return (value: string) => {
    if (type !== TYPES.DOI) return true;

    const isValid = value.includes("doi.org");

    return isValid || "Please enter a valid DOI link";
  };
}

function isValidHttpUrl(url: string): true | string {
  let result;

  try {
    result = new URL(url);
  } catch (_) {
    result = { protocol: "" };
  }

  const isValid = result.protocol === "http:" || result.protocol === "https:";

  return isValid || "Please enter a link";
}

export default LinkInput;
