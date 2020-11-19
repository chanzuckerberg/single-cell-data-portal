import { Button, Classes, Intent } from "@blueprintjs/core";
import { useNavigate } from "@reach/router";
import React, { FC, useEffect, useRef, useState } from "react";
import { ROUTES } from "src/common/constants/routes";
import { COLLECTION_LINK_TYPE } from "src/common/entities";
import {
  formDataToObject,
  useCreateCollection,
} from "src/common/queries/collections";
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
}

type Link = {
  id: number;
  url: string;
  isValid: boolean;
  type: COLLECTION_LINK_TYPE;
};

enum FIELD_NAMES {
  NAME = "name",
  DESCRIPTION = "description",
  CONTACT_NAME = "contact-name",
  CONTACT_EMAIL = "contact-email",
}

const Content: FC<Props> = (props) => {
  const [isValid, setIsValid] = useState(false);
  const [policyVersion, setPolicyVersion] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const navigate = useNavigate();

  const [fieldValidation, setFieldValidation] = useState<{
    [key: string]: boolean;
  }>({
    [FIELD_NAMES.NAME]: false,
    [FIELD_NAMES.DESCRIPTION]: false,
    [FIELD_NAMES.CONTACT_NAME]: false,
    [FIELD_NAMES.CONTACT_EMAIL]: false,
  });

  const formEl = useRef<HTMLFormElement>(null);

  const [links, setLinks] = useState<Link[]>([]);

  const [mutate] = useCreateCollection();

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
      <div className={Classes.DIALOG_BODY}>
        <Form ref={formEl}>
          <Input
            name={FIELD_NAMES.NAME}
            text="Collection Name"
            handleChange={handleInputChange}
            placeholder={REQUIRED_FIELD_TEXT}
          />
          <StyledLabel htmlFor="description">
            <LabelText>Description</LabelText>
            <TextArea
              name={FIELD_NAMES.DESCRIPTION}
              id={FIELD_NAMES.DESCRIPTION}
              handleChange={handleInputChange}
              placeholder={REQUIRED_FIELD_TEXT}
              fill
            />
          </StyledLabel>
          <ContactWrapper>
            <StyledInput
              name={FIELD_NAMES.CONTACT_NAME}
              text="Contact Name"
              handleChange={handleInputChange}
              placeholder={REQUIRED_FIELD_TEXT}
            />
            <StyledInput
              name={FIELD_NAMES.CONTACT_EMAIL}
              text="Contact Email"
              handleChange={handleInputChange}
              placeholder={REQUIRED_FIELD_TEXT}
            />
          </ContactWrapper>
          {links.map(({ type, id }, index) => (
            <LinkInput
              index={index}
              type={type}
              id={id}
              key={id}
              handleChange={handleLinkInputChange}
              handleDelete={handleLinkInputDelete}
            />
          ))}
          <AddLink handleClick={handleAddLinkClick} />
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
            onClick={onSubmit}
            loading={isLoading}
          >
            Create
          </Button>
        </div>
      </div>
    );
  }

  async function onSubmit() {
    if (!formEl?.current) return;

    const formData = new FormData(formEl.current);

    const payload = formDataToObject(formData);

    const payloadLinks = links.map(({ type, url }) => ({
      // eslint-disable-next-line @typescript-eslint/camelcase
      link_type: type,
      // eslint-disable-next-line @typescript-eslint/camelcase
      link_url: url,
    }));

    payload.links = payloadLinks;
    payload[POLICY_PAYLOAD_KEY] = policyVersion;

    setIsLoading(true);

    const collectionId = (await mutate(JSON.stringify(payload))) as string;

    setIsLoading(false);

    if (collectionId) {
      navigate(ROUTES.PRIVATE_COLLECTION.replace(":id", collectionId));
    }
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
    value,
    isValid: isValidFromLinkInput,
  }: LinkValue) {
    const link = links[index];
    const newLink: Link = {
      ...link,
      isValid: isValidFromLinkInput,
      url: value,
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

function createLinkInput(type: COLLECTION_LINK_TYPE) {
  return { id: Date.now(), isValid: false, type, url: "" };
}

export default Content;
