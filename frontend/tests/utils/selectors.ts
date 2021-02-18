export const getTestID = (id: string) => `css=[data-test-id="${id}"]`;

export const getTestClass = (className: string) =>
  `css=[data-test-class="${className}"]`;

export const getText = (text: string) => `text=${text}`;
