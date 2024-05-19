import styled from "@emotion/styled";

export const CardWrapper = styled.div`
  position: fixed;
  bottom: 0;
  left: 50%;
  transform: translateX(-50%);
  width: 90%;
  max-width: 600px;
  margin: 16px;
  max-height: 80vh;
  overflow-y: auto;
  padding: 16px;
  background: #fff;
  box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
  border-radius: 8px;
  z-index: 1000;
  display: flex;
  flex-direction: column;
  align-items: flex-end; /* Align items to the end (right) */
`;

export const CardContent = styled.div`
  font-size: 14px;
  color: #333;
  width: 100%; /* Ensure content takes full width */
`;

export const CloseButtonWrapper = styled.div`
  position: absolute;
  top: 8px;
  right: 8px;
  background-color: white;
`;
