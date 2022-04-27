import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { FC } from "react";
import {
  COLLECTION_LINK_TYPE,
  COLLECTION_LINK_TYPE_OPTIONS,
} from "src/common/entities";
import {
  FormLabelText,
  SelectFormLabel,
} from "src/components/common/Form/common/style";
import Input from "src/components/common/Form/Input";
import { GRAY } from "src/components/common/theme";
import { isLinkTypeDOI } from "src/components/CreateCollectionModal/components/Content/common/utils";
import AddLink from "../AddLink";
import {
  CloseCollectionLinkIcon,
  CollectionLink,
  HelperText,
  InputPrefix,
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
const DOI_PLACEHOLDER = "10.12345/67890123456789";
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
  const option = COLLECTION_LINK_TYPE_OPTIONS[linkType];
  const { text, value } = option;
  const isDOI = isLinkTypeDOI(value);
  const urlPrefix = isDOI ? (
    <InputPrefix warning={!!errorMessage}>doi:</InputPrefix>
  ) : undefined;
  const urlPlaceholder = getUrlPlaceholder(isDOI);

  // Determine validation for link. All links except DOIs are validated on the FE. DOIs are validated by the BE.
  const syncValidation = isDOI ? [] : [isValidHttpUrl];

  const LinkTypeButton = () => (
    <Button fill minimal outlined rightIcon="caret-down" text={text} />
  );

  return (
    <CollectionLink>
      <SelectFormLabel>
        <FormLabelText>Type</FormLabelText>
        <AddLink
          doiSelected={doiSelected}
          fill
          handleClick={handleLinkTypeChange}
          Button={LinkTypeButton}
        />
      </SelectFormLabel>
      {!isDOI && (
        <Input
          name="Name"
          // (thuang): `noNameAttr` removes this input field from the FormData and
          // the payload
          noNameAttr
          optionalField
          text="Name"
          placeholder="Name"
          handleChange={handleNameChange}
          defaultValue={linkName}
        />
      )}
      <Input
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
        isRevalidationRequired={isRevalidationRequired}
      />
      {/* Helper text */}
      {isDOI && (
        <HelperText warning={!!errorMessage}>
          {errorMessage ? errorMessage : DOI_HELPER_TEXT}
        </HelperText>
      )}
      <CloseCollectionLinkIcon
        color={GRAY.A}
        icon={IconNames.CROSS}
        onClick={() => handleDelete(id)}
      />
    </CollectionLink>
  );

  /**
   * Returns "url" input field placeholder text.
   * @param isDoiLink
   * @returns placeholder text for the "url" input field.
   */
  function getUrlPlaceholder(isDoiLink: boolean): string {
    if (isDoiLink) {
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
    // TODO(cc) review comment below; test was isFilterEnabled && isDio (#1718)
    // Adding special handling for filter-enabled DOI: DOI is considered invalid if not specified (that is, prevent form
    // submit) but don't mark the input as an error.
    const isValid = isDOI ? !!value : isValidFromInput;
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
 * @param url - URL specified by user.
 * @param linkType - Currently selected link type.
 * @param newLinkType - Link type to switch to.
 * @param isValid - True if field is currently marked as valid
 * @param isTouched - True if field value has been modified by user.
 * @returns True if re-validation is required.
 */
function isRevalidationRequired_(
  url: string,
  linkType: COLLECTION_LINK_TYPE,
  newLinkType: COLLECTION_LINK_TYPE,
  isValid: boolean,
  isTouched = false
): boolean {
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
