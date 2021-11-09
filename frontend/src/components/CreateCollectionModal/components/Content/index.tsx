import { Button, Classes, Intent } from "@blueprintjs/core";
import { useRouter } from "next/router";
import { FC, useEffect, useRef, useState } from "react";
import { ROUTES } from "src/common/constants/routes";
import {
  Collection,
  COLLECTION_LINK_TYPE,
  VISIBILITY_TYPE,
} from "src/common/entities";
import {
  formDataToObject,
  useCollection,
  useCreateCollection,
  useEditCollection,
} from "src/common/queries/collections";
import { isTombstonedCollection } from "src/common/utils/typeGuards";
import { Value } from "src/components/common/Form/common/constants";
import Input from "src/components/common/Form/Input";
import { LabelText, StyledLabel } from "src/components/common/Form/Input/style";
import TextArea from "src/components/common/Form/TextArea";
import AddLink from "./components/AddLink";
import LinkInput, { LinkValue } from "./components/LinkInput";
import Policy from "./components/Policy";
import { ContactWrapper, Form, StyledInput } from "./style";

const POLICY_PAYLOAD_KEY = "data_submission_policy_version";

const REQUIRED_FIELD_TEXT = "Required";

interface Props {
  onClose: () => void;
  editingMode?: boolean;
  id?: Collection["id"];
}

type Link = {
  id: number;
  url: string;
  linkName: string;
  isValid: boolean;
  linkType: COLLECTION_LINK_TYPE;
};

enum FIELD_NAMES {
  NAME = "name",
  DESCRIPTION = "description",
  CONTACT_NAME = "contact-name",
  CONTACT_EMAIL = "contact-email",
}

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
  const [policyVersion, setPolicyVersion] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const router = useRouter();

  const [fieldValidation, setFieldValidation] = useState<{
    [key: string]: boolean;
  }>({
    [FIELD_NAMES.NAME]: initialBooleanState,
    [FIELD_NAMES.DESCRIPTION]: initialBooleanState,
    [FIELD_NAMES.CONTACT_NAME]: initialBooleanState,
    [FIELD_NAMES.CONTACT_EMAIL]: initialBooleanState,
  });

  const formEl = useRef<HTMLFormElement>(null);

  let { data } = useCollection({
    id: props.id,
    visibility: VISIBILITY_TYPE.PRIVATE,
  });

  if (isTombstonedCollection(data)) data = {} as Collection;

  const [mutateCreateCollection] = useCreateCollection();
  const [mutateEditCollection] = useEditCollection(props.id);

  const { name, description, contact_email, contact_name } = data || {};

  const [links, setLinks] = useState<Link[]>(
    data?.links.map((link, index) => {
      return {
        id: Date.now() + index,
        isValid: true,
        linkName: link.link_name,
        linkType: link.link_type,
        url: link.link_url,
      };
    }) || []
  );

  useEffect(() => {
    const areLinksValid = links.every((link) => link.isValid);
    const areFieldsValid = Object.values(fieldValidation).every(
      (isValid) => isValid
    );
    const isPolicyChecked = policyVersion !== "";

    const result = areLinksValid && areFieldsValid && isPolicyChecked;

    if (result !== isValid) {
      setIsValid(result);
    }
  }, [links, fieldValidation, policyVersion, isValid]);

  const { onClose } = props;
  return (
    <>
      <div
        data-test-id="collection-form-content"
        className={Classes.DIALOG_BODY}
      >
        {/* (thuang): turn off autocomplete seems to be the safest bet for now
         * https://github.com/facebook/react/issues/1159
         */}
        <Form ref={formEl} autoComplete="off">
          <Input
            name={FIELD_NAMES.NAME}
            text="Collection Name"
            handleChange={handleInputChange}
            placeholder={REQUIRED_FIELD_TEXT}
            syncValidation={[requiredValidator]}
            defaultValue={name}
          />
          <StyledLabel htmlFor="description">
            <LabelText>Description</LabelText>
            <TextArea
              name={FIELD_NAMES.DESCRIPTION}
              id={FIELD_NAMES.DESCRIPTION}
              handleChange={handleInputChange}
              placeholder={REQUIRED_FIELD_TEXT}
              fill
              defaultValue={description}
            />
          </StyledLabel>
          <ContactWrapper>
            <StyledInput
              name={FIELD_NAMES.CONTACT_NAME}
              text="Contact Name"
              handleChange={handleInputChange}
              placeholder={REQUIRED_FIELD_TEXT}
              defaultValue={contact_name}
              syncValidation={[requiredValidator]}
            />
            <StyledInput
              name={FIELD_NAMES.CONTACT_EMAIL}
              text="Contact Email"
              handleChange={handleInputChange}
              placeholder={REQUIRED_FIELD_TEXT}
              defaultValue={contact_email}
              syncValidation={[requiredValidator, emailValidation]}
            />
          </ContactWrapper>
          {links.map(({ linkType, id, url, linkName, isValid }, index) => (
            <LinkInput
              index={index}
              linkType={linkType}
              id={id}
              key={id}
              linkName={linkName}
              handleChange={handleLinkInputChange}
              handleDelete={handleLinkInputDelete}
              url={url}
              isValid={isValid}
            />
          ))}
          <AddLink handleClick={handleAddLinkClick} Button={AddLinkButton} />
          <Policy handleChange={handlePolicyChange} />
        </Form>
      </div>
      <Footer />
    </>
  );

  function Footer() {
    return (
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button
            minimal
            intent={Intent.PRIMARY}
            onClick={onClose}
            disabled={isLoading}
          >
            Cancel
          </Button>
          <Button
            intent={Intent.PRIMARY}
            disabled={!isValid}
            onClick={
              isEditCollection ? submitEditCollection : submitCreateCollection
            }
            loading={isLoading}
            data-test-id="create-button"
          >
            {isEditCollection ? "Save" : "Create"}
          </Button>
        </div>
      </div>
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
    payload[POLICY_PAYLOAD_KEY] = policyVersion;

    return payload;
  }

  async function submitCreateCollection() {
    const payload = createPayload();

    if (!payload) return;

    setIsLoading(true);

    const collectionId = (await mutateCreateCollection(
      JSON.stringify(payload)
    )) as string;

    setIsLoading(false);

    if (collectionId) {
      router.push(ROUTES.PRIVATE_COLLECTION.replace(":id", collectionId));
    }
  }

  async function submitEditCollection() {
    const payload = createPayload();

    if (!payload) return;

    setIsLoading(true);

    await mutateEditCollection({
      id: props.id ?? "",
      payload: JSON.stringify(payload),
    });

    setIsLoading(false);
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
    isValid: isValidFromLinkInput,
    linkName,
    linkType,
  }: LinkValue) {
    const link = links[index];
    const newLink: Link = {
      ...link,
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

  function handlePolicyChange(value: string) {
    setPolicyVersion(value);
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
