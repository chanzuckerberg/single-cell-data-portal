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

export const Subtitle = styled.p`
  color: #666;
  margin-bottom: 32px;
  font-size: 16px;
`;

export const FeatureCardsContainer = styled.div`
  display: flex;
  gap: 24px;
  margin-bottom: 58px;
  justify-content: space-between;
  flex-wrap: nowrap;

  @media (max-width: 768px) {
    flex-direction: column;
    flex-wrap: wrap;
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
`;

export const FeatureIconWrapper = styled.div`
  flex-shrink: 0;
  display: flex;
  align-items: flex-start;
  margin-top: 2px;
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
