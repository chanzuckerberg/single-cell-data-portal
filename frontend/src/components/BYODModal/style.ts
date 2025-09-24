import styled from "@emotion/styled";
import { Dialog, DialogContent } from "@czi-sds/components";

export const StyledDialog = styled(Dialog)`
  && .MuiDialog-paper {
    min-height: auto;
    height: auto;
  }

  && .MuiDialogTitle-root {
    padding-right: 64px;

    @media (max-width: 768px) {
      padding-right: 48px;
    }

    @media (max-width: 480px) {
      padding-right: 56px;
      font-size: 1.1rem;
      line-height: 1.3;
    }
  }
`;

export const StyledDialogContent = styled(DialogContent)`
  && {
    height: auto;
    min-height: auto;
    max-height: none;
    overflow: visible;
    flex: 0 0 auto;
    max-width: 1200px;
    width: 100%;
  }
`;

export const FeatureCardsContainer = styled.div`
  display: flex;
  gap: 24px;
  margin-bottom: 58px;
  justify-content: space-between;
  flex-wrap: nowrap;

  /* Tablet: 2 boxes on top, 1 full-width on bottom */
  @media (max-width: 1024px) and (min-width: 769px) {
    flex-wrap: wrap;
    justify-content: flex-start;
  }

  /* Mobile: All boxes stack vertically */
  @media (max-width: 768px) {
    flex-direction: column;
    flex-wrap: nowrap;
  }
`;

export const FeatureCard = styled.div`
  flex: 1;
  min-width: 0;
  border: 1px solid #e0e0e0;
  border-radius: 4px;
  padding: 24px;
  background: white;
  display: flex;
  flex-direction: row;
  gap: 16px;
  box-sizing: border-box;

  /* Tablet: First two cards take ~50% width, third card takes full width */
  @media (max-width: 1024px) and (min-width: 769px) {
    flex: 0 0 calc(50% - 12px);

    &:nth-of-type(3) {
      flex: 1 0 100%;
    }
  }
`;

export const FeatureIconWrapper = styled.div`
  flex-shrink: 0;
  display: flex;
  align-items: flex-start;
`;

export const FeatureContentWrapper = styled.div`
  flex: 1;
  display: flex;
  flex-direction: column;
`;

export const FeatureTitle = styled.h3`
  margin: 0 0 16px 0;
  font-size: 18px;
  font-weight: 600;
  color: #333;
`;

export const FeatureDescription = styled.p`
  color: #666;
  margin: 0 0 24px 0;
  font-size: 14px;
  line-height: 1.5;
  flex: 1;
`;

export const FeatureButtonContainer = styled.div`
  display: flex;
  justify-content: flex-start;
`;

export const ModalFooter = styled.div`
  display: flex;
  justify-content: flex-end;
  align-items: flex-end;
`;

export const StyledIcon = styled.span`
  margin-left: 8px;
  display: inline-flex;
  align-items: center;
`;
