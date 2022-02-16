import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { FC, Fragment } from "react";
import {
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
} from "src/common/entities";
import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import {
  FormLabelText as StyledLabelText,
  SelectFormLabel,
} from "src/components/common/Form/common/style";
import Input from "src/components/common/Form/Input";
import { LabelText, StyledDiv } from "src/components/common/Form/Input/style";
import { GRAY } from "src/components/common/theme";
import AddLink from "../AddLink";
import {
  CloseCollectionLinkIcon,
  CollectionLink,
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
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const LinkInputWrapper = isFilterEnabled ? CollectionLink : LinkWrapper;
  const FormLabel = isFilterEnabled ? SelectFormLabel : StyledDiv;
  const SelectButton = isFilterEnabled ? Button : StyledLinkTypeButton;
  const FormLabelText = isFilterEnabled ? StyledLabelText : LabelText;
  const CloseCollectionLink = isFilterEnabled ? Fragment : IconWrapper;

  // All links except DOIs are validated on the FE. DOIs are validated by the BE.
  // TODO replace syncValidation below with this definition once filter flag is removed.
  const filterEnabledValidation =
    linkType === COLLECTION_LINK_TYPE.DOI ? [] : [isValidHttpUrl];

  // Determine validation for link.
  const syncValidation = isFilterEnabled
    ? filterEnabledValidation
    : [isValidHttpUrl, isDOILink(value)];

  const LinkTypeButton = () => (
    <SelectButton fill minimal outlined rightIcon="caret-down" text={text} />
  );

  return (
    <LinkInputWrapper>
      <FormLabel>
        <FormLabelText>Type</FormLabelText>
        <AddLink
          doiSelected={false}
          fill
          handleClick={handleLinkTypeChange}
          Button={LinkTypeButton}
        />
      </FormLabel>
      <Input
        name="Name"
        // (thuang): `noNameAttr` removes this input field from the FormData and
        // the payload
        noNameAttr
        optionalField={isFilterEnabled}
        text={isFilterEnabled ? "Name" : "Name (optional)"}
        placeholder="Name"
        handleChange={handleNameChange}
        defaultValue={linkName}
        percentage={25} // TODO(cc) remove prop once filter feature flag is removed (#1718).
      />
      <StyledURLInput // TODO(cc) revert to Input component once filter feature flag is removed (#1718).
        // (thuang): `noNameAttr` removes this input field from the FormData and
        // the payload
        noNameAttr
        name={value}
        text="URL"
        syncValidation={syncValidation}
        placeholder={
          value === COLLECTION_LINK_TYPE.DOI
            ? DOI_PLACEHOLDER
            : LINK_PLACEHOLDER
        }
        defaultValue={url}
        handleChange={handleChange_}
        percentage={40} // TODO(cc) remove prop once filter feature flag is removed (#1718).
      />
      <CloseCollectionLink>
        {isFilterEnabled ? (
          <CloseCollectionLinkIcon
            color={GRAY.A}
            icon={IconNames.CROSS}
            onClick={() => handleDelete(id)}
          />
        ) : (
          <StyledButton
            minimal
            color={GRAY.A}
            icon={IconNames.CROSS}
            onClick={() => handleDelete(id)}
          />
        )}
      </CloseCollectionLink>
    </LinkInputWrapper>
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

export function isDOILink(
  type: COLLECTION_LINK_TYPE
): (value: string) => true | string {
  return (value: string) => {
    // Skip validation if type is not DOI
    if (type !== COLLECTION_LINK_TYPE.DOI) return true;
    const origin = "doi.org/";
    const originIndex = value.indexOf("doi.org/");
    const isValid =
      originIndex !== -1 && value.slice(originIndex + origin.length).length > 0;

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
