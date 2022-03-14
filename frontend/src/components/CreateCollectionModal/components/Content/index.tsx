import { Button, Classes, Intent } from "@blueprintjs/core";
import { useRouter } from "next/router";
import { FC, Fragment, useEffect, useRef, useState } from "react";
import { ROUTES } from "src/common/constants/routes";
import {
  Collection,
  COLLECTION_LINK_TYPE,
  VISIBILITY_TYPE,
} from "src/common/entities";
import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { useUserInfo } from "src/common/queries/auth";
import {
  formDataToObject,
  useCollection,
  useCreateCollection,
  useEditCollection,
} from "src/common/queries/collections";
import { isTombstonedCollection } from "src/common/utils/typeGuards";
import {
  StyledDefaultButton,
  StyledPrimaryAnchorButton,
  StyledPrimaryButton,
} from "src/components/common/Button/common/style";
import { Value } from "src/components/common/Form/common/constants";
import Input from "src/components/common/Form/Input";
import TextArea from "src/components/common/Form/TextArea";
import { isLinkTypeDOI } from "src/components/CreateCollectionModal/components/Content/common/utils";
import { getDOIPath } from "src/views/Collection/utils";
import AddLink from "./components/AddLink";
import LinkInput, { LinkValue } from "./components/LinkInput";
import {
  CollectionDetail as Detail,
  CollectionFooter,
  CollectionLinks as Links,
  ContactWrapper,
  Form,
  FormDivider,
  StyledInput,
  Title,
} from "./style";

const REQUIRED_FIELD_TEXT = "Required";

/**
 * Text displayed when BE has identified DOI as invalid.
 */
const INVALID_DOI_ERROR_MESSAGE =
  "This DOI could not be found. Please correct or remove it.";

interface Props {
  onClose: () => void;
  editingMode?: boolean;
  id?: Collection["id"];
}

type Link = {
  errorMessage?: string; // Populated by server-side errors
  id: number;
  url: string;
  linkName: string;
  isRevalidationRequired?: boolean; // True if switching between link fields with different validation (e.g. DOI vs others).
  isTouched?: boolean; // True if field value has been modified by user.
  isValid: boolean;
  linkType: COLLECTION_LINK_TYPE;
};

enum FIELD_NAMES {
  NAME = "name",
  DESCRIPTION = "description",
  CONTACT_NAME = "contact-name",
  CONTACT_EMAIL = "contact-email",
}

const AddMetadataLinkButton = () => (
  <StyledPrimaryAnchorButton intent={Intent.PRIMARY} minimal text="Add Link" />
);

/**
 * @deprecated by AddInputGroupButton once filter feature flag is removed (#1718).
 */
const AddLinkButton = () => (
  <Button outlined intent={Intent.PRIMARY}>
    Add Link
  </Button>
);

const requiredValidator = (value: string) => value.length > 0 || "Required";

// checks if string value is a valid email
const emailValidation = (value: string) => {
  const emailRegex =
    /^(([^<>()[\]\\.,;:\s@"]+(\.[^<>()[\]\\.,;:\s@"]+)*)|(".+"))@((\[[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}])|(([a-zA-Z\-0-9]+\.)+[a-zA-Z]{2,}))$/;

  return emailRegex.test(value) || "Invalid email";
};

const Content: FC<Props> = (props) => {
  const isEditCollection = !!props.id;
  const initialBooleanState = isEditCollection;
  const [isValid, setIsValid] = useState(initialBooleanState);
  const [isLoading, setIsLoading] = useState(false);
  const router = useRouter();
  const { data: userInfo } = useUserInfo(true);
  const [fieldValidation, setFieldValidation] = useState<{
    [key: string]: boolean;
  }>({
    [FIELD_NAMES.NAME]: initialBooleanState,
    [FIELD_NAMES.DESCRIPTION]: initialBooleanState,
    [FIELD_NAMES.CONTACT_NAME]: initialBooleanState,
    [FIELD_NAMES.CONTACT_EMAIL]: initialBooleanState,
  });

  /* Temporary structure to enable both filter feature and existing functionality. */
  /* Fragments and any deprecated styled components can either be removed or replaced, or reverted to bp components once filter feature flag is removed (#1718). */
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const CollectionDetail = isFilterEnabled ? Detail : Fragment;
  const ContactDetailsWrapper = isFilterEnabled
    ? Fragment
    : ContactWrapper; /* remove ContactDetailsWrapper from structure once filter feature flag is removed (#1718) */
  const CollectionLinks = isFilterEnabled ? Links : Fragment;
  const DialogFooter = isFilterEnabled ? CollectionFooter : "div";
  const CancelButton = isFilterEnabled ? StyledDefaultButton : Button;
  const SaveButton = isFilterEnabled ? StyledPrimaryButton : Button;

  const formEl = useRef<HTMLFormElement>(null);

  let { data } = useCollection({
    id: props.id,
    visibility: VISIBILITY_TYPE.PRIVATE,
  });

  const { mutateAsync: mutateCreateCollection } = useCreateCollection();
  const { mutateAsync: mutateEditCollection } = useEditCollection(props.id);

  // Null / tombstone checking is type safety netting.  We shouldn't be getting to these lines/cases since we can't open the modal if the collection is tombstoned/doesn't exist.
  if (isTombstonedCollection(data)) data = null;
  const { name, description, contact_email, contact_name } = data || {};

  const [links, setLinks] = useState<Link[]>([]);
  // TODO useEffect can be reverted back to useState once isFilterEnabled is removed (#1718). See 7210a51d18b2747ea79b47fcad28579bb5b514b6.
  useEffect(() => {
    if (isTombstonedCollection(data)) {
      return;
    }
    setLinks(
      data?.links.map((link, index) => {
        const {
          link_name: linkName,
          link_type: linkType,
          link_url: linkUrl,
        } = link;

        // DOI links are specified as path only, all other links are specified as full URLs.
        let url;
        if (isLinkTypeDOI(linkType)) {
          url = isFilterEnabled ? getDOIPath(linkUrl) : linkUrl;
        } else {
          url = linkUrl;
        }
        return {
          id: Date.now() + index,
          isValid: true,
          linkName,
          linkType,
          url,
        };
      }) || []
    );
  }, [data, isFilterEnabled]);

  useEffect(() => {
    const areLinksValid = links.every((link) => link.isValid);
    const areFieldsValid = Object.values(fieldValidation).every(
      (isValid) => isValid
    );

    const result = areLinksValid && areFieldsValid;

    if (result !== isValid) {
      setIsValid(result);
    }
  }, [links, fieldValidation, isValid]);

  const { onClose } = props;
  return (
    <>
      <div
        data-test-id="collection-form-content"
        className={isFilterEnabled ? undefined : Classes.DIALOG_BODY}
      >
        {/* (thuang): turn off autocomplete seems to be the safest bet for now
         * https://github.com/facebook/react/issues/1159
         */}
        <Form ref={formEl} autoComplete="off" isFilterEnabled={isFilterEnabled}>
          {/* Collection detail */}
          <CollectionDetail>
            {/* Title */}
            {isFilterEnabled && <Title>Collection Details</Title>}
            {/* Fields */}
            <Input
              name={FIELD_NAMES.NAME}
              text="Collection Name"
              handleChange={handleInputChange}
              placeholder={REQUIRED_FIELD_TEXT}
              syncValidation={[requiredValidator]}
              defaultValue={name}
            />
            <TextArea
              name={FIELD_NAMES.DESCRIPTION}
              id={FIELD_NAMES.DESCRIPTION}
              handleChange={handleInputChange}
              placeholder={REQUIRED_FIELD_TEXT}
              fill
              defaultValue={description}
            />
            <ContactDetailsWrapper>
              <StyledInput // @deprecated - revert to bp Input once filter feature flag is removed (#1718).
                name={FIELD_NAMES.CONTACT_NAME}
                text="Contact Name"
                handleChange={handleInputChange}
                placeholder={REQUIRED_FIELD_TEXT}
                defaultValue={contact_name}
                syncValidation={[requiredValidator]}
              />
              <StyledInput // @deprecated - revert to bp Input once filter feature flag is removed (#1718).
                name={FIELD_NAMES.CONTACT_EMAIL}
                text="Contact Email"
                handleChange={handleInputChange}
                placeholder={REQUIRED_FIELD_TEXT}
                defaultValue={contact_email}
                syncValidation={[requiredValidator, emailValidation]}
              />
            </ContactDetailsWrapper>
          </CollectionDetail>
          {/* Collection links */}
          {links.length > 0 && (
            <>
              {isFilterEnabled && <FormDivider />}
              <CollectionLinks>
                {/* Title */}
                {isFilterEnabled && <Title>Links</Title>}
                {/* Fields */}
                {links.map(
                  (
                    {
                      errorMessage,
                      linkType,
                      id,
                      url,
                      linkName,
                      isRevalidationRequired,
                      isTouched,
                      isValid,
                    },
                    index
                  ) => (
                    <LinkInput
                      doiSelected={isDOISelected()}
                      errorMessage={errorMessage}
                      index={index}
                      linkType={linkType}
                      id={id}
                      key={id}
                      linkName={linkName}
                      handleChange={handleLinkInputChange}
                      handleDelete={handleLinkInputDelete}
                      url={url}
                      isRevalidationRequired={isRevalidationRequired}
                      isTouched={isTouched}
                      isValid={isValid}
                    />
                  )
                )}
              </CollectionLinks>
            </>
          )}
          {/* Add metadata link button */}
          <AddLink
            doiSelected={isDOISelected()}
            handleClick={handleAddLinkClick}
            Button={isFilterEnabled ? AddMetadataLinkButton : AddLinkButton}
          />
        </Form>
      </div>
      {/* Modal footer */}
      <Footer />
    </>
  );

  function Footer() {
    return (
      <DialogFooter className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <CancelButton
            disabled={isLoading}
            intent={isFilterEnabled ? undefined : Intent.PRIMARY}
            minimal
            onClick={onClose}
            text="Cancel"
          />
          <SaveButton
            data-test-id="create-button"
            disabled={!isValid}
            intent={Intent.PRIMARY}
            loading={isLoading}
            onClick={
              isEditCollection ? submitEditCollection : submitCreateCollection
            }
            text={isFilterEnabled || isEditCollection ? "Save" : "Create"}
          />
        </div>
      </DialogFooter>
    );
  }

  function createPayload() {
    if (!formEl?.current) return;

    const formData = new FormData(formEl.current);

    const payload = formDataToObject(formData);

    const payloadLinks = links.map(({ linkType, url, linkName: name }) => ({
      link_name: name,
      link_type: linkType,
      link_url: url,
    }));

    payload.links = payloadLinks;
    payload.curator_name = userInfo?.name;
    return payload;
  }

  /**
   * Determine if there is currently a DOI link type selected.
   * @returns True if a DOI link type has been added to the array of links.
   */
  function isDOISelected(): boolean {
    return links.some((link) => isLinkTypeDOI(link.linkType));
  }

  async function submitCreateCollection() {
    const payload = createPayload();

    if (!payload) return;

    setIsLoading(true);

    const { collectionId, isInvalidDOI } = await mutateCreateCollection(
      JSON.stringify(payload)
    );

    setIsLoading(false);

    // Handle the case where DOI is invalid.
    if (isInvalidDOI) {
      setLinks(updateLinkErrors(links, true));
      return;
    }

    if (collectionId) {
      router.push(ROUTES.PRIVATE_COLLECTION.replace(":id", collectionId));
    }
  }

  async function submitEditCollection() {
    const payload = createPayload();

    delete payload?.curator_name; // Do not update curator name when revising a collection

    if (!payload) return;

    setIsLoading(true);

    const response = await mutateEditCollection({
      id: props.id ?? "",
      payload: JSON.stringify(payload),
    });

    setIsLoading(false);

    // Handle the case where DOI is invalid.
    const { isInvalidDOI } = response;
    if (isInvalidDOI) {
      setLinks(updateLinkErrors(links, true));
      return;
    }

    onClose();
  }

  function handleInputChange({ isValid: isValidFromInput, name }: Value) {
    if (fieldValidation[name] === isValidFromInput) return;

    const newFieldValidation = {
      ...fieldValidation,
      [name]: isValidFromInput,
    };

    setFieldValidation(newFieldValidation);
  }

  function handleLinkInputChange({
    index,
    url,
    isRevalidationRequired,
    isTouched,
    isValid: isValidFromLinkInput,
    linkName,
    linkType,
  }: LinkValue) {
    const link = links[index];
    const newLink: Link = {
      ...link,
      errorMessage: "", // Clear server-side errors
      isRevalidationRequired,
      isTouched,
      isValid: isValidFromLinkInput,
      linkName,
      linkType,
      url,
    };

    const newLinks = [...links];
    newLinks[index] = newLink;

    setLinks(newLinks);
  }

  function handleLinkInputDelete(id: number) {
    const newLinks = links.filter((link) => link.id !== id);

    setLinks(newLinks);
  }

  function handleAddLinkClick(type: COLLECTION_LINK_TYPE) {
    const link = createLinkInput(type);

    const newLinks = [...links, link];

    setLinks(newLinks);
  }

  /**
   * Update link validation status to include server-side errors.
   * TODO generalize beyond DOI link type once all links are validated on the BE (#1916).
   * @param links - Current set of selected links.
   * @param isInvalidDOI - True if the server has indicated the submitted DOI is invalid.
   * @returns Array of links with error messages updated according to server-side errors.
   */
  function updateLinkErrors(links: Link[], isInvalidDOI: boolean): Link[] {
    return links.map((link) => {
      if (isLinkTypeDOI(link.linkType) && isInvalidDOI) {
        return {
          ...link,
          errorMessage: INVALID_DOI_ERROR_MESSAGE,
        };
      }
      return link;
    });
  }
};

function createLinkInput(linkType: COLLECTION_LINK_TYPE): Link {
  return {
    id: Date.now(),
    isValid: false,
    linkName: "",
    linkType,
    url: "",
  };
}

export default Content;
