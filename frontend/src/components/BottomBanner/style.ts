import styled from "@emotion/styled";
import {
  Banner,
  Button,
  CommonThemeProps,
  fontBodyS,
  fontBodyXxxs,
  fontHeaderXl,
  getColors,
  InputText,
} from "czifui";
import Modal from "../common/Modal";

export const BOTTOM_BANNER_ID = "bottom-banner";

export const HiddenHubspotForm = styled.div`
  display: none;
`;

export const StyledBanner = styled(Banner)`
  ${fontBodyS}

  padding-left: 16px;

  letter-spacing: -0.006em;

  height: auto;

  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

    // beta intent does not exist for SDS banner, but the colors do
    // targeting specific id to overwrite style
    return `
      border-color: ${colors?.beta[400]} !important;
      background-color: ${colors?.beta[100]};
      color: black;

      /* Hide default svg icon in the Banner as it is not in figma */
      :first-child > div:first-child > div:first-child {
        display: none;
      }

      /* Change close button icon default color */
      button svg {
        path {
          fill: ${colors?.gray[500]};
        }
      }
    `;
  }}
`;

export const StyledBottomBannerWrapper = styled.div`
  ${asFooter}

  width: 100%;

  /* Right behind modal overlay */
  z-index: 19;
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
  width: 100%;
`;

export const StyledTitle = styled.div`
  ${fontHeaderXl}

  width: 100%;

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

  width: 100%;

  letter-spacing: -0.005em;

  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

    // beta intent does not exist for SDS banner, but the colors do
    // targeting specific id to overwrite style
    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;

export const StyledErrorMessage = styled.div`
  ${fontBodyXxxs}

  letter-spacing: -0.005em;

  align-self: flex-start;

  height: 16px;
  margin-top: 4px;
  margin-bottom: 4px;

  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

    // beta intent does not exist for SDS banner, but the colors do
    // targeting specific id to overwrite style
    return `
      color: ${colors?.error[400]};
    `;
  }}
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
    min-height: 100vh !important;
    overflow: hidden !important;
    padding-top: 170px;
    border-radius: 0 !important;
  }
`;

export const FooterContentWrapper = styled.div`
  margin: 24px 24px 40px 24px;
  display: flex;
  flex-direction: column;
  align-items: center;
`;

function asFooter({ asFooter }: { asFooter: boolean }) {
  if (asFooter) {
    return `
      display: block;
    `;
  } else {
    return `
      position: absolute;
      bottom: 0;
    `;
  }
}
