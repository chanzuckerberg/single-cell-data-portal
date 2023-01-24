import { useAuth0 } from "@auth0/auth0-react";
import { AnchorButton, Button, Classes, Intent } from "@blueprintjs/core";
import { FC } from "react";
import { useWindowLocationOrigin } from "src/common/hooks/useWindowLocationOrigin";
import { ContentWrapper } from "./style";

interface Props {
  onClose: () => void;
}

const CTA: FC<Props> = ({ onClose }) => {

  const { loginWithRedirect } = useAuth0();
  const windowLocationOriginUri = useWindowLocationOrigin();

  return (
    <>
      <ContentWrapper>
        Only registered users can create collections.
      </ContentWrapper>
      <Footer />
    </>
  );

  function Footer() {
    return windowLocationOriginUri ? (
      <div className={Classes.DIALOG_FOOTER}>
        <div className={Classes.DIALOG_FOOTER_ACTIONS}>
          <Button minimal intent={Intent.PRIMARY} onClick={onClose}>
            Cancel
          </Button>
          <AnchorButton
            intent={Intent.PRIMARY}
            onClick={() => loginWithRedirect({ redirectUri: `${windowLocationOriginUri}?showCC=true`})}
          >
            Continue
          </AnchorButton>
        </div>
      </div>
    ) : null;
  }
};

export default CTA;
