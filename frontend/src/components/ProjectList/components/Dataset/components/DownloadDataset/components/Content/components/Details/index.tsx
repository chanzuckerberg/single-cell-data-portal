import React, { FC } from "react";
import { Section, Text, Title } from "../common/style";
import { Spinner } from "./style";

interface Props {
  selected: boolean;
  fileSize: number;
  isLoading: boolean;
}

const MEGA_BYTES = 1000 * 1000;

const Details: FC<Props> = ({
  selected = false,
  fileSize = 0,
  isLoading = false,
}) => {
  function renderContent() {
    if (isLoading) {
      return <Spinner size={Spinner.SIZE_SMALL} />;
    }

    if (!selected) {
      return "Select from an above data format to view download details.";
    }

    return `${Math.floor(fileSize / MEGA_BYTES)}MB`;
  }

  return (
    <Section>
      <Title>DOWNLOAD DETAILS</Title>
      <Text>{renderContent()}</Text>
    </Section>
  );
};

export default Details;
