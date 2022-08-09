import { Button, Collapse, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { FC, useState } from "react";
import { GRAY } from "src/components/common/theme";
import DirectIdentifiers from "./components/DirectPersonalIdentifiers";
import {
  BulletWrapper,
  ButtonWrapper,
  ContentWrapper,
  StyledButton,
  StyledIcon,
  Text,
  Wrapper,
} from "./style";

export const POLICY_BULLETS = {
  items: [
    "I give CZI permission to display, distribute, and create derivative works (e.g. visualizations) of this data for purposes of offering CELLxGENE Discover, and I have the authority to give this permission.",
    <div key="1">
      It is my responsibility to ensure that this data is not identifiable. In
      particular, I commit that I will remove any{" "}
      <DirectIdentifiers>
        <ButtonWrapper>
          <StyledButton>direct personal identifiers</StyledButton>
        </ButtonWrapper>
      </DirectIdentifiers>{" "}
      in the metadata portions of the data, and that CZI may further contact me
      if it believes more work is needed to de-identify it.
    </div>,
    <div key="2">
      If I choose to publish this data <b>publicly</b> on CELLxGENE Discover, I
      understand that (1) anyone will be able to access it subject to a CC-BY
      license, meaning they can download, share, and use the data without
      restriction beyond providing attribution to the original data
      contributor(s) and (2) the Collection details (including collection name,
      description, my name, and the contact information for the datasets in this
      Collection) will be made public on the CELLxGENE Discover as well.
    </div>,
    "I understand that I have the ability to delete the data that I have published from CELLxGENE Discover if I later choose to. This however will not undo any prior downloads or shares of such data.",
  ],
  version: "2.0",
};

const Policy: FC = () => {
  const [isOpen, setIsOpen] = useState(false);

  return (
    <Wrapper>
      By publishing this collection, you agree to CELLxGENE&#39;s data
      submission policies.
      <Button minimal intent={Intent.PRIMARY} onClick={handleShowButtonClick}>
        {isOpen ? "Hide" : "Show"} Details
      </Button>
      <Collapse isOpen={isOpen}>
        <ContentWrapper>
          {POLICY_BULLETS.items.map((item, index) => (
            <Bullet key={index} content={item} />
          ))}
        </ContentWrapper>
      </Collapse>
    </Wrapper>
  );

  function handleShowButtonClick() {
    setIsOpen(!isOpen);
  }
};

function Bullet({ content }: { content: string | JSX.Element }) {
  return (
    <BulletWrapper>
      <StyledIcon color={GRAY.A} iconSize={10} icon={IconNames.SYMBOL_CIRCLE} />
      <Text>{content}</Text>
    </BulletWrapper>
  );
}

export default Policy;
