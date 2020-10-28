import React, { FC } from "react";
import { Section, Text, Title } from "../common/style";

interface Props {
  name: string;
}

const Name: FC<Props> = ({ name }) => (
  <Section>
    <Title>NAME</Title>
    <Text data-test-id="download-asset-name">{name}</Text>
  </Section>
);

export default Name;
