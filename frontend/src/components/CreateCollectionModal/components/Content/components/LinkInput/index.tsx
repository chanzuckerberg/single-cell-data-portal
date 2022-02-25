import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { FC, Fragment } from "react";
import {
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
  COLLECTION_LINK_TYPE_OPTIONS_DEPRECATED,
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
import { isLinkTypeDOI } from "src/components/CreateCollectionModal/components/Content/common/utils";
import AddLink from "../AddLink";
import {
  CloseCollectionLinkIcon,
  CollectionLink,
  HelperText,
  IconWrapper,
  InputPrefix,
  LinkWrapper,
  StyledButton,
  StyledLinkTypeButton,
  StyledURLInput,
} from "./style";

export type LinkValue = {
  id: number;
  url: string;
  isRevalidationRequired?: boolean; // True if switching between link fields with different validation (e.g. DOI vs others).
  isTouched?: boolean; // True if field value has been modified by user.
  isValid: boolean;
  index: number;
  linkName: string;
  linkType: COLLECTION_LINK_TYPE;
};

const DOI_HELPER_TEXT =
  "A summary citation linked to this DOI will be automatically added to this collection.";
const DOI_PLACEHOLDER = "https://doi.org/10.1126/science.aax6234";
const FILTER_ENABLED_DOI_PLACEHOLDER = "10.12345/67890123456789";
const LINK_PLACEHOLDER = "https://cellxgene.cziscience.com";

interface Props {
  doiSelected: boolean;
  errorMessage?: string; // Populated from server side errors
  handleChange: ({ id, url, isValid }: LinkValue) => void;
  id: number;
  index: number;
  isRevalidationRequired?: boolean;
  isTouched?: boolean;
  linkName: string;
  linkType: COLLECTION_LINK_TYPE;
  handleDelete: (id: number) => void;
  url: string;
  isValid: boolean;
}

const LinkInput: FC<Props> = ({
  doiSelected,
  errorMessage,
  handleChange,
  handleDelete,
  id,
  isRevalidationRequired,
  linkType,
  index,
  url,
  linkName,
  isTouched,
  isValid,
}) => {
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const option = isFilterEnabled
    ? COLLECTION_LINK_TYPE_OPTIONS[linkType]
    : COLLECTION_LINK_TYPE_OPTIONS_DEPRECATED[linkType];
  const { text, value } = option;
  const LinkInputWrapper = isFilterEnabled ? CollectionLink : LinkWrapper;
  const FormLabel = isFilterEnabled ? SelectFormLabel : StyledDiv;
  const SelectButton = isFilterEnabled ? Button : StyledLinkTypeButton;
  const FormLabelText = isFilterEnabled ? StyledLabelText : LabelText;
  const CloseCollectionLink = isFilterEnabled ? Fragment : IconWrapper;
  const isDOI = isLinkTypeDOI(value);
  const isFilterEnabledDOI = isFilterEnabled && isDOI;
  const nameFieldVisible = isNameFieldVisible(isFilterEnabled, isDOI);
  const urlPrefix = isFilterEnabledDOI ? (
    <InputPrefix warning={!!errorMessage}>doi:</InputPrefix>
  ) : undefined;
  const urlPlaceholder = getUrlPlaceholder(isFilterEnabled, isDOI);

  // All links except DOIs are validated on the FE. DOIs are validated by the BE.
  // TODO replace syncValidation below with this definition once filter flag is removed.
  const filterEnabledValidation = isDOI ? [] : [isValidHttpUrl];

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
          doiSelected={doiSelected}
          fill
          handleClick={handleLinkTypeChange}
          Button={LinkTypeButton}
        />
      </FormLabel>
      {nameFieldVisible && (
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
          percentage={25} // TODO remove prop once filter feature flag is removed (#1718).
        />
      )}
      <StyledURLInput // TODO revert to Input component once filter feature flag is removed (#1718).
        // (thuang): `noNameAttr` removes this input field from the FormData and
        // the payload
        leftElement={urlPrefix}
        markAsError={!!errorMessage}
        noNameAttr
        name={value}
        text="URL"
        syncValidation={syncValidation}
        placeholder={urlPlaceholder}
        defaultValue={url}
        handleChange={handleChange_}
        percentage={40} // TODO remove prop once filter feature flag is removed (#1718).
        isRevalidationRequired={isRevalidationRequired}
      />
      {/* Helper text */}
      {isFilterEnabledDOI && (
        <HelperText warning={!!errorMessage}>
          {errorMessage ? errorMessage : DOI_HELPER_TEXT}
        </HelperText>
      )}
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

  /**
   * Returns "url" input field placeholder text.
   * @param isFilterEnabled
   * @param isDoiLink
   * @returns placeholder text for the "url" input field.
   */
  function getUrlPlaceholder(
    isFilterEnabled: boolean,
    isDoiLink: boolean
  ): string {
    if (isDoiLink) {
      if (isFilterEnabled) {
        return FILTER_ENABLED_DOI_PLACEHOLDER;
      }
      return DOI_PLACEHOLDER;
    }
    return LINK_PLACEHOLDER;
  }

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
    // Check if revalidation of link is required.
    const isRevalidationRequired = isRevalidationRequired_(
      isFilterEnabled,
      url,
      linkType,
      newLinkType,
      isValid,
      isTouched
    );
    handleChange({
      id,
      index,
      isRevalidationRequired,
      isTouched,
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
    // Adding special handling for filter-enabled DOI: DOI is considered invalid if not specified (that is, prevent form
    // submit) but don't mark the input as an error.
    const isValid = isFilterEnabledDOI ? !!value : isValidFromInput;
    handleChange({
      id,
      index,
      isTouched: true,
      isValid,
      linkName,
      linkType,
      url: value,
    });
  }

  /**
   * Returns true if either:
   * - filter feature flag is not enabled, or
   * - filter feature flag is enabled and link type is not DOI.
   * @param isFilterEnabled
   * @param isDoiLink
   * @returns boolean
   */
  function isNameFieldVisible(isFilterEnabled: boolean, isDoiLink: boolean) {
    if (!isFilterEnabled) {
      return true;
    }
    return !isDoiLink;
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

/**
 * Determine if validation of field needs to be executed on change on link type.  This is only true if:
 * 1. Filter functionality is enabled.
 * 2. The current link type is DOI or the new link type is DOI (all other link types share the same validation
 * and revalidation is therefore not required).
 * 3. There is currently a value specified for the link URL, or the field has been modified by the user and is
 * currently marked as invalid.
 * @param isFilterEnabled - True if filter flag is enabled.
 * @param url - URL specified by user.
 * @param linkType - Currently selected link type.
 * @param newLinkType - Link type to switch to.
 * @param isValid - True if field is currently marked as valid
 * @param isTouched - True if field value has been modified by user.
 * @returns True if re-validation is required.
 */
function isRevalidationRequired_(
  isFilterEnabled: boolean,
  url: string,
  linkType: COLLECTION_LINK_TYPE,
  newLinkType: COLLECTION_LINK_TYPE,
  isValid: boolean,
  isTouched = false
): boolean {
  // Revalidation is only required for filter-related functionality.
  if (!isFilterEnabled) {
    return false;
  }

  // If neither the current link type nor new link type is DOI, revalidation is not required.
  const isCurrentLinkTypeDOI = isLinkTypeDOI(linkType);
  if (!isCurrentLinkTypeDOI && !isLinkTypeDOI(newLinkType)) {
    return false;
  }

  // Revalidation is not required if current link type is a DOI and there is no value specified.
  if (isCurrentLinkTypeDOI && !url) {
    return false;
  }

  // Revalidate if there is a URL or if the field has been modified and is invalid.
  return !!url || (isTouched && !isValid);
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
