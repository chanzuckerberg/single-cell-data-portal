import styled from "@emotion/styled";
import {
  Banner,
  Button,
  ButtonIcon,
  fontBodyS,
  fontBodyXxxs,
  fontHeaderXl,
  InputText,
} from "@czi-sds/components";
import Modal from "../common/Modal";
import { beta100, beta400, error400, gray500 } from "src/common/theme";
import { SKINNY_MODE_BREAKPOINT_WIDTH } from "src/views/CellGuide/components/CellGuideCard/constants";

export const BOTTOM_BANNER_ID = "bottom-banner";

export const HiddenHubSpotForm = styled.div`
  display: none;
`;

export const StyledBanner = styled(Banner)`
  ${fontBodyS}

  letter-spacing: -0.006em;

  height: auto;

  @media (max-width: ${SKINNY_MODE_BREAKPOINT_WIDTH}px) {
    padding: 8px 16px;
    box-shadow: 0px 0px 4px 0px rgba(50, 50, 50, 0.75);
  }

  /**
   * beta intent does not exist for SDS banner, but the colors do targeting
   * specific id to overwrite style
   */
  border-color: ${beta400} !important;
  background-color: ${beta100};
  color: black;

  /* Hide default svg icon in the Banner as it is not in figma */
  :first-child > div:first-child > div:first-child {
    display: none;
  }

  /* Change close button icon default color */
  button svg {
    path {
      fill: ${gray500};
    }
  }
`;

export const StyledBottomBannerWrapper = styled.div`
  ${asFooter}

  width: 100%;

  /* Right behind modal overlay */
  z-index: 19;

  background-color: purple;
`;

export const StyledLink = styled.a`
  text-decoration-line: underline;
  color: #8f5aff;
  font-weight: 500;

  :hover {
    color: #5826c1;
  }
`;

export const HeaderContainer = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  align-items: flex-start;
`;

const STYLED_CLOSE_BUTTON_ICON_DENY_PROPS = ["hideCloseButton"];

export const StyledCloseButtonIcon = styled(ButtonIcon, {
  shouldForwardProp: (prop) =>
    !STYLED_CLOSE_BUTTON_ICON_DENY_PROPS.includes(prop),
})`
  /* Only hide close button if mobile view and is a direct link */
  @media only screen and (max-width: 600px) {
    ${hideCloseButton}
  }
`;

export const StyledTitle = styled.div`
  ${fontHeaderXl}

  letter-spacing: -0.019em;
  font-size: 24px !important;
  margin: 0;
  height: auto !important;

  padding-top: 16px;
  padding-bottom: 8px;
`;

export const StyledDescription = styled.div`
  ${fontBodyS}

  letter-spacing: -0.006em;
  padding-bottom: 16px;
`;

export const StyledForm = styled.form`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  height: 34px;
  margin-bottom: 0;
  align-items: center;
  width: 100%;
`;

export const StyledInputText = styled(InputText)`
  .MuiOutlinedInput-root.Mui-focused .MuiOutlinedInput-notchedOutline {
    border-color: #5826c1 !important;
  }

  flex: 1;
  margin-right: 4px;
  margin-bottom: 0px;
  display: inline-flex;
`;

export const StyledSubmitButton = styled(Button)`
  padding: 6px 12px;
  width: 91px;
  height: 34px;
  background: #8f5aff;
  font-weight: 500;

  :hover {
    background: #5826c1;
  }
`;

export const StyledDisclaimer = styled.div`
  ${fontBodyXxxs}

  letter-spacing: -0.005em;

  /*
   * beta intent does not exist for SDS banner, but the colors do
   * targeting specific id to overwrite style
   */
  color: ${gray500};
`;

export const StyledErrorMessage = styled.div`
  ${fontBodyXxxs}

  letter-spacing: -0.005em;

  align-self: flex-start;

  height: 16px;
  margin-top: 4px;
  margin-bottom: 4px;

  /*
   * beta intent does not exist for SDS banner, but the colors do
   * targeting specific id to overwrite style
   */
  color: ${error400};
`;

export const NewsletterModal = styled(Modal)`
  .bp4-dialog-header {
    display: none !important;
  }

  min-width: 400px !important;
  min-height: 266px !important;
  max-width: 400px !important;
  max-height: 266px !important;

  margin: 0;

  padding: 24px;

  padding-bottom: 24px !important;

  @media only screen and (max-width: 600px) {
    min-width: 100vw !important;
    max-width: 100vw !important;
    min-height: 100vh !important;
    max-height: 100vh !important;
    overflow: hidden !important;
    padding-top: 170px;
    border-radius: 0 !important;
    position: fixed;
    top: 0;
    left: 0;

    /* Hack to disable touch scrolling when mobile modal is opened */
    touch-action: none;
  }
`;

export const FooterContentWrapper = styled.div`
  margin: 24px 24px 40px 24px;
  display: flex;
  flex-direction: column;
  align-items: center;
`;

function asFooter({ asFooter }: { asFooter: boolean }) {
  // If rendered as part of the footer then have it as a block element that goes with the flow of the page, instead of having it "floating"
  if (asFooter) {
    return `
      display: block;
      text-align: center;
    `;
  } else {
    return `
      position: fixed;
      bottom: 0;
    `;
  }
}

function hideCloseButton({ hideCloseButton }: { hideCloseButton: boolean }) {
  if (!hideCloseButton) return null;
  return `
    visibility: hidden;
  `;
}
