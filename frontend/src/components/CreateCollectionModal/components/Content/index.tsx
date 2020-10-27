import { Button, Classes, Intent } from "@blueprintjs/core";
import { debounce } from "lodash-es";
import React, { FC, useRef, useState } from "react";
import { useMutation, useQueryCache } from "react-query";
import {
  createCollection,
  formDataToObject,
  USE_COLLECTIONS,
} from "src/common/queries/collections";
import Input from "./components/Input";
import { LabelText, StyledLabel } from "./components/Input/style";
import { ContactWrapper, Form, StyledInputGroup } from "./style";

const DEBOUNCE_TIME_MS = 200;

interface Props {
  onClose: () => void;
}

const Content: FC<Props> = (props) => {
  const [isValid, setIsValid] = useState(false);

  const formEl = useRef<HTMLFormElement>(null);

  const nameEl = useRef<HTMLInputElement>(null);
  const descriptionEl = useRef<HTMLInputElement>(null);
  const contactNameEl = useRef<HTMLInputElement>(null);
  const contactEmailEl = useRef<HTMLInputElement>(null);

  const queryCache = useQueryCache();

  const [mutate] = useMutation(createCollection, {
    onSuccess: () => {
      queryCache.invalidateQueries(USE_COLLECTIONS);
    },
  });

  const debouncedValidateFields = debounce(validateFields, DEBOUNCE_TIME_MS);

  const { onClose } = props;

  return (
    <>
      <div className={Classes.DIALOG_BODY}>
        <Form ref={formEl}>
          <Input
            inputRef={nameEl}
            name="name"
            text="Collection Name"
            onChange={debouncedValidateFields}
          />
          <Input
            inputRef={descriptionEl}
            name="description"
            text="Description"
            onChange={debouncedValidateFields}
          />
          <ContactWrapper>
            <StyledLabel htmlFor="contact-name">
              <LabelText>Contact Name</LabelText>
              <StyledInputGroup
                inputRef={contactNameEl}
                name="contact-name"
                id="contact-name"
                onChange={debouncedValidateFields}
              />
            </StyledLabel>
            <Input
              inputRef={contactEmailEl}
              name="contact-email"
              text="Contact Email"
              onChange={debouncedValidateFields}
            />
          </ContactWrapper>
        </Form>
      </div>
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button minimal intent={Intent.PRIMARY} onClick={onClose}>
            Cancel
          </Button>
          <Button
            intent={Intent.PRIMARY}
            disabled={!isValid}
            onClick={onSubmit}
          >
            Create
          </Button>
        </div>
      </div>
    </>
  );

  function onSubmit() {
    // payload add links
    // mutate(JSON.string(payload))
    // set collectionId and redirect
    // link_type: DOI, RAW_DATA, PROTOCOL, LAB_WEBSITE, OTHER

    if (!formEl?.current) return;

    const formData = new FormData(formEl.current);

    const payload = formDataToObject(formData);

    // DEBUG
    // DEBUG
    // DEBUG
    console.log("---payload", payload);
  }

  function validateFields() {
    if (
      nameEl.current?.value.length &&
      descriptionEl.current?.value.length &&
      contactNameEl.current?.value.length &&
      contactEmailEl.current?.value.length
    ) {
      return setIsValid(true);
    }

    if (isValid) setIsValid(false);
  }
};

export default Content;
