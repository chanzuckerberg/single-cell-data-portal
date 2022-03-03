import { Dialog, Intent, MenuItem } from "@blueprintjs/core";
import {
  Button as SDSButton,
  DialogActions,
  DialogContent,
  DialogTitle,
} from "czifui";
import { useCallback, useRef, useState } from "react";
import {
  APIKeyResponse,
  generateCuratorAuthKey,
  useCuratorAuthKeyID,
} from "src/common/queries/curation";
import Toast from "src/views/Collection/components/Toast";
import { APIDisclaimerP, FullWidthCallout, StyledInputGroup } from "./style";

export default function CuratorAPIKeyGenerator(): JSX.Element {
  const { data: apiKeyID } = useCuratorAuthKeyID();

  const apiKeyResponse = useRef({ id: apiKeyID ?? "" } as APIKeyResponse);

  const [isOpen, setIsOpen] = useState(false);
  const handleOpen = useCallback(async () => {
    try {
      const data = await generateCuratorAuthKey();
      if (!data?.key) throw new Error("No key returned");
      apiKeyResponse.current = data;
      setIsOpen(!isOpen);
    } catch (error) {
      Toast.show({
        intent: Intent.DANGER,
        message: "An error occurred while generating an API Key",
      });
    }
  }, [isOpen]);
  const handleClose = useCallback(() => setIsOpen(false), []);

  const copyAPIKey = useCallback(() => {
    navigator.clipboard.writeText(apiKeyResponse.current.key || "");
    Toast.show({
      intent: Intent.SUCCESS,
      message: "API Key copied to clipboard",
    });
    setIsOpen(false);
  }, [apiKeyResponse]);

  return (
    <>
      <MenuItem
        data-testid="get-api-key"
        text={`${
          apiKeyResponse.current.id !== "" ? "Regenerate" : "Generate"
        } API Key`}
        onClick={handleOpen}
        shouldDismissPopover={false}
      />
      <Dialog
        isOpen={isOpen}
        onClose={handleClose}
        canOutsideClickClose={false}
      >
        <DialogTitle title="API Key" />
        <DialogContent>
          <FullWidthCallout intent="info">
            Make sure to copy your API token now. You won’t be able to see it
            again!
          </FullWidthCallout>
          <StyledInputGroup readOnly value={apiKeyResponse.current.key} />
          <APIDisclaimerP>
            Warning: Treat your tokens like passwords and keep them secret. When
            working with the API, use tokens as environment variables instead of
            hardcoding them into your programs.
          </APIDisclaimerP>
        </DialogContent>
        <DialogActions>
          <SDSButton onClick={handleClose}>Close</SDSButton>
          <SDSButton sdsType="primary" sdsStyle="square" onClick={copyAPIKey}>
            Copy API Key
          </SDSButton>
        </DialogActions>
      </Dialog>
    </>
  );
}
