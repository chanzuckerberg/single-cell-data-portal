import { Intent } from "@blueprintjs/core";
/* eslint-disable @typescript-eslint/ban-ts-comment --- TODO revisit lint errors below */
import {
  Button as SDSButton,
  // @ts-ignore
  Dialog,
  // @ts-ignore
  DialogActions,
  // @ts-ignore
  DialogContent,
  // @ts-ignore
  DialogTitle,
  MenuItem,
} from "czifui";
/* eslint-enable @typescript-eslint/ban-ts-comment --- TODO revisit lint errors above */
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
  const handleClose = useCallback((_, reason) => {
    if (reason === "backdropClick") return;

    setIsOpen(false);
  }, []);

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
      <MenuItem data-testid="get-api-key" onClick={handleOpen}>
        {apiKeyResponse.current.id !== "" ? "Regenerate" : "Generate"} API Key
      </MenuItem>
      <Dialog open={isOpen} onClose={handleClose}>
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
          {/* eslint-disable-next-line @typescript-eslint/ban-ts-comment -- TODO revisit lint error */}
          {/* @ts-ignore -- TODO revisit lint error */}
          <SDSButton onClick={handleClose}>Close</SDSButton>
          <SDSButton sdsType="primary" sdsStyle="square" onClick={copyAPIKey}>
            Copy API Key
          </SDSButton>
        </DialogActions>
      </Dialog>
    </>
  );
}
