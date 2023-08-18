import { Intent } from "@blueprintjs/core";
import {
  Button as SDSButton,
  Dialog,
  DialogActions,
  DialogContent,
  DialogProps,
  DialogTitle,
  MenuItem,
} from "@czi-sds/components";
import { useCallback, useRef, useState } from "react";
import {
  APIKeyResponse,
  generateCuratorAuthKey,
  useCuratorAuthKeyID,
} from "src/common/queries/curation";
import Toast from "src/views/Collection/components/Toast";
import { APIDisclaimerP, FullWidthCallout, StyledInputGroup } from "./style";

type DialogOnCloseParams = Parameters<NonNullable<DialogProps["onClose"]>>;

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

  const handleClose = useCallback(
    (_: DialogOnCloseParams[0], reason: DialogOnCloseParams[1]) => {
      if (reason === "backdropClick") return;

      setIsOpen(false);
    },
    []
  );

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
            Make sure to copy your API key now. You won’t be able to see it
            again!
          </FullWidthCallout>
          <StyledInputGroup readOnly value={apiKeyResponse.current.key} />
          <APIDisclaimerP>
            Warning: Treat your API key like a password and keep it secret.
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
